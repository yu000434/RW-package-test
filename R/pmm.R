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
  model_names <- model$xnames
  predictors <- setdiff(model_names, "(Intercept)")
  x <- as.matrix(data[predictors])
  if ("(Intercept)" %in% model_names) x <- cbind(`(Intercept)` = 1, x)
  x <- x[, model_names, drop = FALSE]
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

#' Compute the PMM cross term
#'
#' @param object An `rw_fit` object returned by [with_rw()].
#' @param variable Name of the PMM-imputed variable.
#' @param expected_score A function called with `object`, `variable`, the
#'   imputation number, the missing-row indicator, and the donor predictions.
#'   It must return one donor-by-recipient matrix per analysis coefficient.
#' @param ... Additional arguments passed to `expected_score`.
#'
#' @return The PMM columns of the RW cross-term matrix.
#' @export
pmm_kappa <- function(object, variable, expected_score, ...) {
  imps <- object$mids
  analysis_names <- names(stats::coef(object$results[[1L]]$model))
  pmm_names <- imps$models[[variable]][[1L]]$xnames
  kappa <- matrix(0, length(analysis_names), length(pmm_names),
                  dimnames = list(analysis_names, pmm_names))

  for (p in seq_len(object$m)) {
    inputs <- pmm_inputs(imps, variable, p)
    score <- expected_score(
      object, variable, p, inputs$missing,
      inputs$membership$donor_prediction, ...
    )
    expected_dim <- dim(inputs$membership$probability)
    if (length(score) != length(analysis_names) ||
        any(vapply(score, function(x) !identical(dim(x), expected_dim), logical(1)))) {
      stop("`expected_score` returned incompatible score matrices.")
    }
    donors <- imps$models[[variable]][[p]]$setup$donors
    for (j in seq_along(pmm_names)) {
      derivative <- inputs$membership$derivative[[j]] / donors
      kappa[, j] <- kappa[, j] +
        vapply(score, function(x) sum(derivative * x), numeric(1))
    }
  }
  kappa / (object$n * object$m)
}

pmm_score_lm <- function(object, variable, p, missing, donor_mean) {
  model <- object$results[[p]]$model
  if (!inherits(model, "lm") || inherits(model, "glm") ||
      all.vars(stats::formula(model))[[1L]] != variable) {
    stop("Automatic PMM pooling requires the PMM variable to be the response of lm().")
  }
  coef_names <- names(stats::coef(model))
  x <- matrix(0, object$n, length(coef_names), dimnames = list(NULL, coef_names))
  rows <- as.integer(rownames(stats::model.frame(model)))
  x[rows, ] <- stats::model.matrix(model)
  x <- x[missing, , drop = FALSE]
  residual <- outer(donor_mean, drop(x %*% stats::coef(model)), "-")
  lapply(seq_along(coef_names), function(j) {
    sweep(residual, 2L, x[, j] / summary(model)$sigma^2, "*")
  })
}

pmm_lm_kappa <- function(object, variable) {
  pmm_kappa(object, variable, pmm_score_lm)
}

gauss_legendre <- function(order) {
  index <- seq_len(order - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  eig <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(eig$values)
  list(node = eig$values[ordering], weight = 2 * eig$vectors[1L, ordering]^2)
}

pmm_score_binomial <- function(object, variable, p, missing, donor_mean,
                               threshold, quadrature_order) {
  model <- object$results[[p]]$model
  response <- all.vars(stats::formula(model))[[1L]]
  analysis_names <- names(stats::coef(model))
  if (!inherits(model, "glm") || stats::family(model)$family != "binomial" ||
      !identical(analysis_names, c("(Intercept)", variable))) {
    stop("Binomial PMM pooling requires glm(response ~ PMM_variable, family = binomial()).")
  }
  completed <- mice::complete(object$mids, p)
  included <- as.integer(rownames(stats::model.frame(model)))
  if (!identical(included, which(completed[[variable]] > threshold))) {
    stop("The analysis subset must be the PMM variable greater than `threshold`.")
  }

  downstream <- object$mids$models[[response]][[p]]
  if (is.null(downstream) || object$mids$method[[response]] != "logreg") {
    stop("The binomial response must use logreg imputation.")
  }
  beta <- drop(downstream$beta.dot)
  model_names <- downstream$xnames
  model_names[model_names == ""] <- "(Intercept)"
  names(beta) <- model_names
  base_names <- setdiff(model_names, variable)
  predictors <- setdiff(base_names, "(Intercept)")
  x <- as.matrix(object$mids$data[missing, predictors, drop = FALSE])
  if ("(Intercept)" %in% base_names) x <- cbind(`(Intercept)` = 1, x)
  downstream_base <- drop(x[, base_names, drop = FALSE] %*% beta[base_names])
  observed_response <- object$mids$data[[response]][missing]
  if (is.factor(observed_response)) {
    observed_response <- as.numeric(as.character(observed_response))
  }

  donor_sd <- object$mids$models[[variable]][[p]]$sigma.dot
  rule <- gauss_legendre(quadrature_order)
  lower <- pmax((threshold - donor_mean) / donor_sd, -10)
  upper <- rep(10, length(donor_mean))
  half_width <- ifelse(lower < upper, (upper - lower) / 2, 0)
  midpoint <- (upper + lower) / 2
  z <- outer(half_width, rule$node, "*") + midpoint
  weight <- outer(half_width, rule$weight, "*") * stats::dnorm(z)
  donor_value <- sweep(z * donor_sd, 1L, donor_mean, "+")
  score <- list(matrix(0, length(donor_mean), sum(missing)),
                matrix(0, length(donor_mean), sum(missing)))
  observed <- !is.na(observed_response)

  for (q in seq_len(quadrature_order)) {
    value <- donor_value[, q]
    analysis_mean <- stats::plogis(stats::coef(model)[[1L]] +
                                     stats::coef(model)[[variable]] * value)
    downstream_mean <- stats::plogis(outer(value * beta[[variable]], downstream_base, "+"))
    residual <- sweep(downstream_mean, 1L, analysis_mean, "-")
    if (any(observed)) {
      residual[, observed] <- outer(analysis_mean, observed_response[observed],
                                    function(mean, outcome) outcome - mean)
    }
    weighted <- residual * weight[, q]
    score[[1L]] <- score[[1L]] + weighted
    score[[2L]] <- score[[2L]] + sweep(weighted, 1L, value, "*")
  }
  score
}

#' Compute the PMM cross term for a binomial analysis
#'
#' @param object An `rw_fit` object returned by [with_rw()].
#' @param variable Name of the continuous PMM predictor.
#' @param threshold Lower analysis threshold for the PMM predictor.
#' @param quadrature_order Number of quadrature points.
#'
#' @return The PMM columns of the RW cross-term matrix.
#' @export
pmm_kappa_binomial <- function(object, variable, threshold = -Inf,
                               quadrature_order = 24L) {
  pmm_kappa(object, variable, pmm_score_binomial,
            threshold = threshold, quadrature_order = quadrature_order)
}
