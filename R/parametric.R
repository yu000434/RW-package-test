parametric_component <- function(data, model, variable) {
  beta <- setNames(as.numeric(model$beta.dot), model$xnames)
  names(beta)[names(beta) == ""] <- "(Intercept)"
  terms <- setdiff(names(beta), "(Intercept)")
  x <- as.matrix(data[terms])
  if ("(Intercept)" %in% names(beta)) x <- cbind(`(Intercept)` = 1, x)
  x <- x[, names(beta), drop = FALSE]
  y <- data[[variable]]
  imputed <- data[[paste0(".imputed_", variable)]]
  observed <- !imputed

  if (model$setup$method == "logreg") {
    mean <- plogis(drop(x %*% beta))
    score <- x * (y - mean)
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
    information <- matrix(0, p + 1L, p + 1L)
    information[seq_len(p), seq_len(p)] <- -crossprod(x_obs) / sigma2
    information[seq_len(p), p + 1L] <- -drop(crossprod(x_obs, residual_obs)) / sigma2^2
    information[p + 1L, seq_len(p)] <- information[seq_len(p), p + 1L]
    information[p + 1L, p + 1L] <-
      sum(1 / (2 * sigma2^2) - residual_obs^2 / sigma2^3)
    information <- information / nrow(data)
  }

  list(S_mis_imp = score * imputed,
       d = t(-solve(information, t(score * observed))))
}
