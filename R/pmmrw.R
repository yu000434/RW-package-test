#' Predictive mean matching with donor recording
#'
#' Performs standard MICE predictive mean matching while recording donor
#' and model information needed for Robins-Wang variance estimation.
#'
#' @param y,ry,x,wy Arguments supplied by `mice`.
#' @param task,model Model-recording arguments supplied by `mice`.
#' @param ... Additional arguments.
#'
#' @return Imputed values from the selected donors.
#' @export
#' 
mice.impute.pmmrw <- function(y, ry, x, wy = !ry, task = "impute", model = NULL,
                              exclude = NULL, ridge = 1e-05, matchtype = 1L,
                              donors = 5L, use.matcher = FALSE, mlocal = 1L, ...) {
  if (!is.null(exclude) || isTRUE(use.matcher) || mlocal != 1L || matchtype != 1L) {
    stop("`pmmrw` requires exclude = NULL, use.matcher = FALSE, mlocal = 1, and matchtype = 1.")
  }
  if (is.factor(y)) stop("`pmmrw` supports numeric variables only.")
  if (task != "train") stop("`pmmrw` requires `tasks = \"train\"`.")

  n <- length(y)
  x <- cbind(`(Intercept)` = 1, as.matrix(x))
  draw <- utils::getFromNamespace(".norm.draw", "mice")(y, ry, x, ridge = ridge, ...)
  beta_hat <- drop(draw$coef)
  beta_dot <- drop(draw$beta)

  x_obs <- x[ry, , drop = FALSE]
  x_mis <- x[wy, , drop = FALSE]
  pred_obs <- drop(x_obs %*% beta_hat)
  pred_mis <- drop(x_mis %*% beta_dot)
  if (is.null(donors)) donors <- round(length(pred_obs) / 600 + 7)
  donors <- max(1L, min(as.integer(donors), length(pred_obs)))
  donor_pos <- mice::matchindex(pred_obs, pred_mis, donors)

  sigma <- as.numeric(draw$sigma)
  residual <- y[ry] - drop(x_obs %*% draw$coef)
  information <- -crossprod(x_obs) / (sigma^2 * n)
  d_obs <- t(-solve(information, t(x_obs * (residual / sigma^2))))
  pmm_d <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_d[ry, ] <- d_obs
  donor_id <- rep(NA_integer_, n)
  donor_id[wy] <- which(ry)[donor_pos]

  model$setup <- list(method = "pmmrw", n = sum(ry), task = task, donors = donors,
                      matchtype = matchtype, ridge = ridge)
  model$beta.hat <- beta_hat
  model$beta.dot <- beta_dot
  model$sigma.dot <- sigma
  model$xnames <- colnames(x)
  model$donor_id <- donor_id
  model$pmm_score <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  model$pmm_d <- pmm_d

  y[ry][donor_pos]
}
