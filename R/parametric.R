parametric_component <- function(data, model, variable, imputed) {
  # Extract the imputation-model parameter draw beta and align coefficient names.
  beta <- setNames(as.numeric(model$beta.dot), model$xnames)
  names(beta)[names(beta) == ""] <- "(Intercept)"
  
  # Design matrix for the imputation model.
  x <- imputation_design(data, model)
  # Completed outcome: contains observed values and imputed values.
  y <- data[[variable]]
  # Indicator for whether each value of y was originally missing and then imputed.
  observed <- !imputed

  if (model$setup$method == "logreg") {
    if (is.factor(y)) y <- as.integer(y) - 1L
    mean <- plogis(drop(x %*% beta))
    # Individual score contribution for the imputation model.
    score <- x * (y - mean)
    # Score Jacobian, using only originally observed outcomes.
    information <- -crossprod(x[observed, , drop = FALSE] *
                                sqrt(mean[observed] * (1 - mean[observed]))) / nrow(data)
  } else {
    sigma2 <- model$sigma.dot^2
    residual <- y - drop(x %*% beta)
    score <- cbind(x * (residual / sigma2),
                   sigma2 = 0.5 * (-1 / sigma2 + residual^2 / sigma2^2))
    x_obs <- x[observed, , drop = FALSE]
    residual_obs <- residual[observed]
    p <- ncol(x)
    # Jacobian of the score with respect to (beta, sigma^2).
    information <- matrix(0, p + 1L, p + 1L)
    information[seq_len(p), seq_len(p)] <- -crossprod(x_obs) / sigma2
    information[seq_len(p), p + 1L] <- -drop(crossprod(x_obs, residual_obs)) / sigma2^2
    information[p + 1L, seq_len(p)] <- information[seq_len(p), p + 1L]
    information[p + 1L, p + 1L] <-
      sum(1 / (2 * sigma2^2) - residual_obs^2 / sigma2^3)
    information <- information / nrow(data)
  }

  list(
    # Score contributions from originally missing/imputed observations.
    S_mis_imp = score * imputed,
    # Influence contribution
    d = t(-solve(information, t(score * observed))))
}
