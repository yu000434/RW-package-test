# Constructs analysis components for linear and generalized linear models.

analysis_component <- function(model, n) {
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  residual <- y - fitted(model)

  if (inherits(model, "glm")) {
    U <- x * residual
    mu <- predict(model, type = "response")
    link <- predict(model, type = "link")
    weight <- family(model)$mu.eta(link)^2 / family(model)$variance(mu)
    tau <- -crossprod(x * sqrt(weight))
  } else {
    sigma2 <- summary(model)$sigma^2
    U <- x * residual / sigma2
    tau <- -crossprod(x) / sigma2
  }

  full_U <- matrix(0, n, ncol(U), dimnames = list(NULL, names(coef(model))))
  full_U[as.integer(rownames(model.frame(model))), ] <- U
  list(U = full_U, tau = tau, n_analysis = nrow(U))
}
