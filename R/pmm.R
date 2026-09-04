pmm_component <- function(model) {
  list(S_mis_imp = model$pmm_score, d = model$pmm_d)
}

#' Extract recorded PMM donor IDs
#'
#' @param data A `mids` object containing a `pmmrw` imputation.
#' @param variable Name of the PMM-imputed variable.
#'
#' @return An `n` by `m` matrix of donor row IDs.
#' @export
extract_donor_id <- function(data, variable) {
  do.call(cbind, lapply(seq_len(data$m), function(p) {
    data$models[[variable]][[p]]$donor_id
  }))
}

topk_membership <- function(x_obs, x_mis, beta_hat, beta_dot, covariance, donors) {
  x_obs <- as.matrix(x_obs)
  x_mis <- as.matrix(x_mis)
  n_obs <- nrow(x_obs)
  n_par <- ncol(x_obs)
  if (donors < 1L || donors > n_obs) stop("Invalid donor pool size.")

  donor_prediction <- drop(x_obs %*% beta_hat)
  target <- drop(x_mis %*% beta_dot)
  smooth_sd <- sqrt(rowSums((x_mis %*% covariance) * x_mis))
  if (any(!is.finite(smooth_sd) | smooth_sd <= 0)) stop("Invalid matching variance.")
  ordering <- order(donor_prediction)
  sorted_prediction <- donor_prediction[ordering]
  sorted_x <- x_obs[ordering, , drop = FALSE]
  rank <- seq_len(n_obs)
  has_lower <- rank > donors
  has_upper <- rank + donors <= n_obs

  lower_sorted <- rep(-Inf, n_obs)
  upper_sorted <- rep(Inf, n_obs)
  lower_slope_sorted <- upper_slope_sorted <- matrix(0, n_obs, n_par)
  lower_sorted[has_lower] <- (sorted_prediction[has_lower] +
    sorted_prediction[rank[has_lower] - donors]) / 2
  upper_sorted[has_upper] <- (sorted_prediction[has_upper] +
    sorted_prediction[rank[has_upper] + donors]) / 2
  lower_slope_sorted[has_lower, ] <- (sorted_x[has_lower, , drop = FALSE] +
    sorted_x[rank[has_lower] - donors, , drop = FALSE]) / 2
  upper_slope_sorted[has_upper, ] <- (sorted_x[has_upper, , drop = FALSE] +
    sorted_x[rank[has_upper] + donors, , drop = FALSE]) / 2

  lower <- upper <- numeric(n_obs)
  lower_slope <- upper_slope <- matrix(0, n_obs, n_par)
  lower[ordering] <- lower_sorted
  upper[ordering] <- upper_sorted
  lower_slope[ordering, ] <- lower_slope_sorted
  upper_slope[ordering, ] <- upper_slope_sorted

  z_lower <- sweep(outer(lower, target, "-"), 2L, smooth_sd, "/")
  z_upper <- sweep(outer(upper, target, "-"), 2L, smooth_sd, "/")
  probability <- stats::pnorm(z_upper) - stats::pnorm(z_lower)
  density_lower <- stats::dnorm(z_lower)
  density_upper <- stats::dnorm(z_upper)
  derivative <- lapply(seq_len(n_par), function(j) {
    upper_difference <- outer(upper_slope[, j], x_mis[, j], "-")
    lower_difference <- outer(lower_slope[, j], x_mis[, j], "-")
    sweep(density_upper * upper_difference - density_lower * lower_difference,
          2L, smooth_sd, "/")
  })
  names(derivative) <- colnames(x_obs)

  tolerance <- 1e-8 * max(1, donors)
  probability_error <- max(abs(colSums(probability) - donors))
  derivative_error <- max(vapply(derivative, function(x) max(abs(colSums(x))), numeric(1)))
  if (probability_error > tolerance || derivative_error > tolerance) {
    stop("Top-k matching probabilities failed their normalization check.")
  }
  intercept <- match("(Intercept)", colnames(x_obs))
  intercept_error <- if (is.na(intercept)) NA_real_ else max(abs(derivative[[intercept]]))
  if (!is.na(intercept_error) && intercept_error > tolerance) {
    stop("A common intercept shift changed a matching probability.")
  }

  list(probability = probability, derivative = derivative,
       donor_prediction = donor_prediction,
       probability_error = probability_error,
       derivative_error = derivative_error,
       intercept_error = intercept_error)
}

pmm_inputs <- function(imps, variable, p) {
  data <- imps$data
  model <- imps$models[[variable]][[p]]
  if (!identical(as.integer(model$setup$matchtype), 1L)) {
    stop("PMM variance estimation requires matchtype = 1.")
  }
  observed <- !is.na(data[[variable]])
  names <- model$xnames
  predictors <- setdiff(names, "(Intercept)")
  x <- as.matrix(data[predictors])
  if ("(Intercept)" %in% names) x <- cbind(`(Intercept)` = 1, x)
  x <- x[, names, drop = FALSE]
  x_obs <- x[observed, , drop = FALSE]
  x_mis <- x[!observed, , drop = FALSE]
  fit <- utils::getFromNamespace("estimice", "mice")(
    x_obs, data[[variable]][observed], ridge = model$setup$ridge
  )
  membership <- topk_membership(
    x_obs, x_mis, drop(model$beta.hat), drop(model$beta.dot),
    model$sigma.dot^2 * fit$v, model$setup$donors
  )
  list(missing = !observed, membership = membership)
}

pmm_columns <- function(imps, variable) {
  variables <- names(imps$models)
  variables <- variables[vapply(variables, function(x) anyNA(imps$data[[x]]), logical(1))]
  dimensions <- vapply(variables, function(x) {
    model <- imps$models[[x]][[1L]]
    switch(unname(imps$method[[x]]),
           pmmrw = ncol(model$pmm_d),
           norm = length(model$beta.dot) + 1L,
           logreg = length(model$beta.dot))
  }, integer(1))
  start <- cumsum(c(1L, dimensions[-length(dimensions)]))
  index <- match(variable, variables)
  seq.int(start[index], length.out = dimensions[index])
}

pmm_lm_kappa <- function(object, variable) {
  n <- object$n
  imps <- object$mids
  analysis_names <- names(stats::coef(object$results[[1L]]$model))
  pmm_names <- imps$models[[variable]][[1L]]$xnames
  kappa <- matrix(0, length(analysis_names), length(pmm_names),
                  dimnames = list(analysis_names, pmm_names))

  for (p in seq_len(object$m)) {
    model <- object$results[[p]]$model
    if (!inherits(model, "lm") || inherits(model, "glm") ||
        all.vars(stats::formula(model))[[1L]] != variable) {
      stop("Automatic PMM pooling requires the PMM variable to be the response of lm().")
    }
    inputs <- pmm_inputs(imps, variable, p)
    x <- matrix(0, n, length(analysis_names), dimnames = list(NULL, analysis_names))
    rows <- as.integer(rownames(stats::model.frame(model)))
    x[rows, ] <- stats::model.matrix(model)
    x_mis <- x[inputs$missing, , drop = FALSE]
    fitted_mean <- drop(x_mis %*% stats::coef(model))
    residual <- outer(inputs$membership$donor_prediction, fitted_mean, "-")
    for (j in seq_along(pmm_names)) {
      weighted <- inputs$membership$derivative[[j]] * residual /
        imps$models[[variable]][[p]]$setup$donors
      kappa[, j] <- kappa[, j] +
        drop(crossprod(x_mis, colSums(weighted))) / summary(model)$sigma^2
    }
  }
  kappa / (n * object$m)
}
