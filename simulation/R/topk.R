# Constructs soft top-k donor probabilities and their PMM-coefficient
# derivatives. These probabilities affect variance estimation only.

topk_membership <- function(
    x_obs,
    x_mis,
    beta_hat,
    beta_dot,
    covariance,
    donors) {
  x_obs <- as.matrix(x_obs)
  x_mis <- as.matrix(x_mis)
  n_obs <- nrow(x_obs)
  n_par <- ncol(x_obs)
  if (ncol(x_mis) != n_par) {
    stop("Observed and missing design matrices have incompatible columns.")
  }
  if (donors < 1L || donors > n_obs) {
    stop("`donors` must be between 1 and the number of observed donors.")
  }

  donor_prediction <- drop(x_obs %*% beta_hat)
  target <- drop(x_mis %*% beta_dot)
  smooth_sd <- sqrt(rowSums((x_mis %*% covariance) * x_mis))
  if (any(!is.finite(smooth_sd) | smooth_sd <= 0)) {
    stop("Every smoothing standard deviation must be positive and finite.")
  }

  ordering <- order(donor_prediction)
  sorted_prediction <- donor_prediction[ordering]
  sorted_x <- x_obs[ordering, , drop = FALSE]
  ranks <- seq_len(n_obs)
  has_lower <- ranks > donors
  has_upper <- ranks + donors <= n_obs

  lower_sorted <- rep(-Inf, n_obs)
  upper_sorted <- rep(Inf, n_obs)
  lower_slope_sorted <- matrix(0, n_obs, n_par)
  upper_slope_sorted <- matrix(0, n_obs, n_par)
  # Each donor's interval is determined by the midpoint to donors k ranks away.
  lower_sorted[has_lower] <- (
    sorted_prediction[has_lower] +
      sorted_prediction[ranks[has_lower] - donors]
  ) / 2
  upper_sorted[has_upper] <- (
    sorted_prediction[has_upper] +
      sorted_prediction[ranks[has_upper] + donors]
  ) / 2
  lower_slope_sorted[has_lower, ] <- (
    sorted_x[has_lower, , drop = FALSE] +
      sorted_x[ranks[has_lower] - donors, , drop = FALSE]
  ) / 2
  upper_slope_sorted[has_upper, ] <- (
    sorted_x[has_upper, , drop = FALSE] +
      sorted_x[ranks[has_upper] + donors, , drop = FALSE]
  ) / 2

  lower <- upper <- numeric(n_obs)
  lower_slope <- upper_slope <- matrix(0, n_obs, n_par)
  lower[ordering] <- lower_sorted
  upper[ordering] <- upper_sorted
  lower_slope[ordering, ] <- lower_slope_sorted
  upper_slope[ordering, ] <- upper_slope_sorted

  z_lower <- sweep(outer(lower, target, "-"), 2L, smooth_sd, "/")
  z_upper <- sweep(outer(upper, target, "-"), 2L, smooth_sd, "/")
  density_lower <- dnorm(z_lower)
  density_upper <- dnorm(z_upper)
  probability <- pnorm(z_upper) - pnorm(z_lower)
  # Move donor predictions, interval boundaries, and recipient predictions
  # together when differentiating with respect to the PMM coefficients.
  derivative <- lapply(seq_len(n_par), function(j) {
    upper_difference <- outer(upper_slope[, j], x_mis[, j], "-")
    lower_difference <- outer(lower_slope[, j], x_mis[, j], "-")
    sweep(
      density_upper * upper_difference -
        density_lower * lower_difference,
      2L,
      smooth_sd,
      "/"
    )
  })
  names(derivative) <- colnames(x_obs)

  probability_error <- max(abs(colSums(probability) - donors))
  derivative_error <- max(vapply(
    derivative,
    function(x) max(abs(colSums(x))),
    numeric(1)
  ))
  tolerance <- 1e-8 * max(1, donors)
  if (probability_error > tolerance || derivative_error > tolerance) {
    stop("Top-k membership failed its normalization check.")
  }

  intercept <- match("(Intercept)", colnames(x_obs))
  intercept_error <- if (is.na(intercept)) {
    NA_real_
  } else {
    max(abs(derivative[[intercept]]))
  }
  if (!is.na(intercept_error) && intercept_error > tolerance) {
    stop("A common intercept shift changed a matching probability.")
  }

  list(
    probability = probability,
    derivative = derivative,
    donor_prediction = donor_prediction,
    target = target,
    smooth_sd = smooth_sd,
    lower = lower,
    upper = upper,
    lower_slope = lower_slope,
    upper_slope = upper_slope,
    probability_error = probability_error,
    derivative_error = derivative_error,
    intercept_error = intercept_error
  )
}

pmm_inputs <- function(
    imps,
    data,
    variable,
    predictors,
    donors,
    p) {
  observed <- !is.na(data[[variable]])
  missing <- !observed
  x <- cbind(`(Intercept)` = 1, as.matrix(data[predictors]))
  x_obs <- x[observed, , drop = FALSE]
  x_mis <- x[missing, , drop = FALSE]
  model <- imps$models[[variable]][[p]]
  if (!identical(as.integer(model$setup$matchtype), 1L)) {
    stop("The joint-boundary prototype currently requires matchtype = 1.")
  }

  estimice <- utils::getFromNamespace("estimice", "mice")
  fit <- estimice(
    x_obs,
    data[[variable]][observed],
    ridge = model$setup$ridge
  )
  covariance <- model$sigma.dot^2 * fit$v
  membership <- topk_membership(
    x_obs,
    x_mis,
    drop(model$beta.hat),
    drop(model$beta.dot),
    covariance,
    donors
  )

  list(
    observed = observed,
    missing = missing,
    x_obs = x_obs,
    x_mis = x_mis,
    y_obs = data[[variable]][observed],
    covariance = covariance,
    membership = membership
  )
}

pmm_cols <- function(imps, variable) {
  imputed_vars <- names(imps$models)
  imputed_vars <- imputed_vars[vapply(
    imputed_vars,
    function(x) anyNA(imps$data[[x]]),
    logical(1)
  )]
  dimensions <- vapply(imputed_vars, function(x) {
    method <- unname(imps$method[[x]])
    model <- imps$models[[x]][[1]]
    if (method == "pmmrw") {
      ncol(model$pmm_d)
    } else if (method == "norm") {
      length(model$beta.dot) + 1L
    } else if (method == "logreg") {
      length(model$beta.dot)
    } else {
      stop("Unsupported imputation method in prototype: ", method)
    }
  }, integer(1))
  start <- cumsum(c(1L, head(dimensions, -1L)))
  index <- match(variable, imputed_vars)
  if (is.na(index)) {
    stop("The requested PMM variable was not imputed.")
  }
  seq.int(start[index], length.out = dimensions[index])
}
