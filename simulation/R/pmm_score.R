# Constructs PMM matching probabilities and their coefficient derivatives.

topk_membership <- function(x_obs, x_mis, beta_hat, beta_dot, covariance, donors) {
  x_obs <- as.matrix(x_obs)
  x_mis <- as.matrix(x_mis)
  n_obs <- nrow(x_obs)
  n_par <- ncol(x_obs)
  if (ncol(x_mis) != n_par) stop("Observed and missing design matrices have incompatible columns.")
  if (donors < 1L || donors > n_obs) stop("`donors` must be between 1 and observed donors.")

  donor_prediction <- drop(x_obs %*% beta_hat)
  target <- drop(x_mis %*% beta_dot)
  smooth_sd <- sqrt(rowSums((x_mis %*% covariance) * x_mis))
  if (any(!is.finite(smooth_sd) | smooth_sd <= 0)) stop("Invalid smoothing standard deviation.")

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
  lower_sorted[has_lower] <- (sorted_prediction[has_lower] +
    sorted_prediction[ranks[has_lower] - donors]) / 2
  upper_sorted[has_upper] <- (sorted_prediction[has_upper] +
    sorted_prediction[ranks[has_upper] + donors]) / 2
  lower_slope_sorted[has_lower, ] <- (sorted_x[has_lower, , drop = FALSE] +
    sorted_x[ranks[has_lower] - donors, , drop = FALSE]) / 2
  upper_slope_sorted[has_upper, ] <- (sorted_x[has_upper, , drop = FALSE] +
    sorted_x[ranks[has_upper] + donors, , drop = FALSE]) / 2

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
    sweep(density_upper * upper_difference - density_lower * lower_difference,
          2L, smooth_sd, "/")
  })
  names(derivative) <- colnames(x_obs)

  probability_error <- max(abs(colSums(probability) - donors))
  derivative_error <- max(vapply(derivative, function(x) max(abs(colSums(x))), numeric(1)))
  tolerance <- 1e-8 * max(1, donors)
  if (probability_error > tolerance || derivative_error > tolerance) {
    stop("Top-k membership failed its normalization check.")
  }

  intercept <- match("(Intercept)", colnames(x_obs))
  intercept_error <- if (is.na(intercept)) NA_real_ else max(abs(derivative[[intercept]]))
  if (!is.na(intercept_error) && intercept_error > tolerance) {
    stop("A common intercept shift changed a matching probability.")
  }

  list(
    probability = probability,
    derivative = derivative,
    donor_prediction = donor_prediction,
    probability_error = probability_error,
    derivative_error = derivative_error,
    intercept_error = intercept_error
  )
}

pmm_inputs <- function(imps, data, variable, predictors, donors, p) {
  observed <- !is.na(data[[variable]])
  missing <- !observed
  x <- cbind(`(Intercept)` = 1, as.matrix(data[predictors]))
  x_obs <- x[observed, , drop = FALSE]
  x_mis <- x[missing, , drop = FALSE]
  model <- imps$models[[variable]][[p]]
  if (!identical(as.integer(model$setup$matchtype), 1L)) stop("PMM score requires matchtype = 1.")

  estimice <- utils::getFromNamespace("estimice", "mice")
  fit <- estimice(x_obs, data[[variable]][observed], ridge = model$setup$ridge)
  covariance <- model$sigma.dot^2 * fit$v
  membership <- topk_membership(x_obs, x_mis, drop(model$beta.hat), drop(model$beta.dot),
                                covariance, donors)

  list(missing = missing, x_mis = x_mis, membership = membership)
}

pmm_cols <- function(imps, variable) {
  imputed_vars <- names(imps$models)
  imputed_vars <- imputed_vars[vapply(imputed_vars, function(x) anyNA(imps$data[[x]]), logical(1))]
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
      stop("Unsupported imputation method: ", method)
    }
  }, integer(1))
  start <- cumsum(c(1L, head(dimensions, -1L)))
  index <- match(variable, imputed_vars)
  if (is.na(index)) stop("The requested PMM variable was not imputed.")
  seq.int(start[index], length.out = dimensions[index])
}

# Computes the PMM cross term for the RW and GS analyses.

pmm_kappa_rw <- function(imps, fit, data, variable, predictors, donors, scenario) {
  n <- nrow(data)
  kappa_sum <- matrix(0, 1L, length(predictors) + 1L,
                      dimnames = list("X", c("(Intercept)", predictors)))
  probability_error <- derivative_error <- intercept_error <- numeric(imps$m)

  for (p in seq_len(imps$m)) {
    inputs <- pmm_inputs(imps, data, variable, predictors, donors, p)
    analysis_model <- fit$results[[p]]$model
    analysis_x <- data$X[inputs$missing]
    analysis_beta <- unname(coef(analysis_model)[["X"]])
    analysis_sigma2 <- summary(analysis_model)$sigma^2

    # E(U_ij | PMM predictors) for each possible donor-recipient pair.
    conditional_residual <- outer(inputs$membership$donor_prediction,
                                  analysis_x * analysis_beta, "-")
    conditional_score <- sweep(conditional_residual, 2L, analysis_x / analysis_sigma2, "*")
    if (scenario == "robins_1") {
      excluded <- data$A[inputs$missing] != 1
      conditional_score[, excluded] <- 0
    }

    kappa_p <- matrix(0, 1L, ncol(inputs$x_mis),
                      dimnames = list("X", colnames(inputs$x_mis)))
    for (component in seq_len(ncol(inputs$x_mis))) {
      derivative <- inputs$membership$derivative[[component]] / donors
      kappa_p[1L, component] <- sum(derivative * conditional_score)
    }
    kappa_sum <- kappa_sum + kappa_p
    probability_error[[p]] <- inputs$membership$probability_error
    derivative_error[[p]] <- inputs$membership$derivative_error
    intercept_error[[p]] <- inputs$membership$intercept_error
  }

  diagnostics <- c(
    max_membership_probability_error = max(probability_error),
    max_membership_derivative_error = max(derivative_error),
    max_intercept_derivative = max(intercept_error)
  )
  list(kappa = kappa_sum / (n * imps$m), diagnostics = diagnostics)
}

# The GS score requires one-dimensional normal integration over donor values.
gauss_legendre_rule <- function(order) {
  index <- seq_len(order - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  eig <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(eig$values)
  list(node = eig$values[ordering], weight = 2 * eig$vectors[1L, ordering]^2)
}

gs_expected_score <- function(donor_mean, donor_sd, downstream_base, downstream_slope,
                              analysis_coefficient, threshold, observed_downstream,
                              quadrature_order = 32L, normal_bound = 10) {
  rule <- gauss_legendre_rule(quadrature_order)
  lower <- pmax((threshold - donor_mean) / donor_sd, -normal_bound)
  upper <- rep(normal_bound, length(donor_mean))
  active <- lower < upper
  half_width <- ifelse(active, (upper - lower) / 2, 0)
  midpoint <- (upper + lower) / 2
  z <- outer(half_width, rule$node, "*") + midpoint
  integration_weight <- outer(half_width, rule$weight, "*") * dnorm(z)
  donor_value <- sweep(z * donor_sd, 1L, donor_mean, "+")

  expected_intercept_score <- matrix(0, length(donor_mean), length(downstream_base))
  expected_slope_score <- expected_intercept_score
  observed <- !is.na(observed_downstream)
  for (node_index in seq_len(quadrature_order)) {
    value <- donor_value[, node_index]
    analysis_mean <- plogis(analysis_coefficient[[1L]] + analysis_coefficient[[2L]] * value)
    downstream_mean <- plogis(outer(value * downstream_slope, downstream_base, "+"))
    residual <- sweep(downstream_mean, 1L, analysis_mean, "-")
    if (any(observed)) {
      residual[, observed] <- outer(analysis_mean, observed_downstream[observed],
                                    function(mean, outcome) outcome - mean)
    }
    weighted_residual <- residual * integration_weight[, node_index]
    expected_intercept_score <- expected_intercept_score + weighted_residual
    expected_slope_score <- expected_slope_score +
      sweep(weighted_residual, 1L, value, "*")
  }

  list(score = list(expected_intercept_score, expected_slope_score),
       omitted_tail_bound = 2 * pnorm(-normal_bound))
}

pmm_kappa_gs <- function(imps, fit, data, variable, predictors, donors, threshold = 2,
                         quadrature_order = 32L, normal_bound = 10) {
  n <- nrow(data)
  analysis_names <- names(coef(fit$results[[1L]]$model))
  pmm_names <- c("(Intercept)", predictors)
  kappa_sum <- matrix(0, length(analysis_names), length(pmm_names),
                      dimnames = list(analysis_names, pmm_names))
  probability_error <- derivative_error <- intercept_error <- numeric(imps$m)
  tail_bound <- numeric(imps$m)

  for (p in seq_len(imps$m)) {
    inputs <- pmm_inputs(imps, data, variable, predictors, donors, p)
    analysis_model <- fit$results[[p]]$model
    analysis_coefficient <- coef(analysis_model)

    downstream <- imps$models$D[[p]]
    downstream_beta <- drop(downstream$beta.dot)
    downstream_names <- downstream$xnames
    downstream_names[downstream_names == ""] <- "(Intercept)"
    names(downstream_beta) <- downstream_names
    base_names <- setdiff(downstream_names, variable)
    downstream_x <- cbind(`(Intercept)` = 1,
                          as.matrix(data[inputs$missing, predictors, drop = FALSE]))
    downstream_base <- drop(downstream_x[, base_names, drop = FALSE] %*%
                              downstream_beta[base_names])
    observed_downstream <- data$D[inputs$missing]
    if (is.factor(observed_downstream)) {
      observed_downstream <- as.numeric(as.character(observed_downstream))
    }

    # Integrate the downstream logistic score under the fitted PMM model.
    conditional <- gs_expected_score(
      inputs$membership$donor_prediction,
      imps$models[[variable]][[p]]$sigma.dot,
      downstream_base,
      downstream_beta[[variable]],
      analysis_coefficient,
      threshold,
      observed_downstream,
      quadrature_order,
      normal_bound
    )

    kappa_p <- matrix(0, length(analysis_names), ncol(inputs$x_mis),
                      dimnames = list(analysis_names, colnames(inputs$x_mis)))
    for (component in seq_len(ncol(inputs$x_mis))) {
      derivative <- inputs$membership$derivative[[component]] / donors
      kappa_p[, component] <- vapply(conditional$score,
                                     function(score) sum(derivative * score), numeric(1))
    }
    kappa_sum <- kappa_sum + kappa_p

    probability_error[[p]] <- inputs$membership$probability_error
    derivative_error[[p]] <- inputs$membership$derivative_error
    intercept_error[[p]] <- inputs$membership$intercept_error
    tail_bound[[p]] <- conditional$omitted_tail_bound
  }

  diagnostics <- c(
    max_membership_probability_error = max(probability_error),
    max_membership_derivative_error = max(derivative_error),
    max_intercept_derivative = max(intercept_error),
    max_quadrature_tail_bound = max(tail_bound)
  )
  list(kappa = kappa_sum / (n * imps$m), diagnostics = diagnostics)
}
