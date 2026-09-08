#' @importFrom stats coef df.residual family fitted model.frame model.matrix
#'   model.response plogis pnorm predict pt qnorm qt setNames
NULL


# Extract the estimating equation components from a fitted analysis model.

analysis_component <- function(model, n) {
  # Design matrix X and outcome vector Y for observations used in fitting the model.
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  # Convert a binary factor outcome to 0/1 for binomial GLMs.
  if (inherits(model, "glm") && is.factor(y)) y <- as.integer(y) - 1L
  residual <- y - fitted(model)

  if (inherits(model, "glm")) {
    if (family(model)$family != "binomial") stop("Only binomial glm models are supported.")
    
    # estimating-function contribution for the regression coefficient beta.
    U <- x * residual
    mu <- predict(model, type = "response")
    link <- predict(model, type = "link")
    weight <- family(model)$mu.eta(link)^2 / family(model)$variance(mu)
    tau <- -crossprod(x * sqrt(weight))
    
  } else if (inherits(model, "lm")) {
    # Estimated residual variance in the Gaussian linear model.
    sigma2 <- summary(model)$sigma^2
    U <- x * residual / sigma2
    tau <- -crossprod(x) / sigma2
  } else {
    stop("The analysis must return an lm or binomial glm model.")
  }
  
  # Expand U back to an n x p matrix aligned with the original dataset.
  #
  # Observations excluded from the fitted model receive a zero score contribution.
  full_U <- matrix(0, n, ncol(U), dimnames = list(NULL, names(coef(model))))
  full_U[as.integer(rownames(model.frame(model))), ] <- U
  list(U = full_U, tau = tau, n_analysis = nrow(U))
}
