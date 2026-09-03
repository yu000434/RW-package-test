# Defines `pmmrw`, a PMM implementation that calls `mice::matchindex()` 
# and records the realized donor IDs and PMM model quantities needed for variance estimation.

mice.impute.pmmrw <- function(
    y,
    ry,
    x,
    wy = NULL,
    task = "impute",
    model = NULL,
    exclude = NULL,
    ridge = 1e-05,
    matchtype = 1L,
    donors = 5L,
    use.matcher = FALSE,
    mlocal = 1L,
    ...) {
  if (!is.null(exclude)) {
    stop("`pmmrw` does not currently support `exclude`.", call. = FALSE)
  }
  if (is.factor(y)) {
    stop("`pmmrw` currently supports numeric variables only.", call. = FALSE)
  }
  if (isTRUE(use.matcher)) {
    stop("`pmmrw` currently supports `use.matcher = FALSE` only.", call. = FALSE)
  }
  if (!identical(as.integer(mlocal), 1L)) {
    stop("`pmmrw` currently supports `mlocal = 1` only.", call. = FALSE)
  }
  if (is.null(wy)) {
    wy <- !ry
  }
  if (task != "train") {
    stop("`pmmrw` must be used with `tasks = 'train'`.", call. = FALSE)
  }
  if (is.null(model) || !is.environment(model)) {
    stop("`model` must be an environment; use `tasks = 'train'`.", call. = FALSE)
  }

  n <- length(y)
  row_id <- seq_len(n)
  x <- cbind(`(Intercept)` = 1, as.matrix(x))
  norm_draw <- utils::getFromNamespace(".norm.draw", "mice")
  parm <- norm_draw(y, ry, x, ridge = ridge, ...)
  beta_hat <- drop(parm$coef)
  beta_dot <- drop(parm$beta)

  if (matchtype == 0L) {
    beta_dot <- beta_hat
  } else if (matchtype == 2L) {
    beta_hat <- beta_dot
  }

  x_obs <- x[ry, , drop = FALSE]
  x_mis <- x[wy, , drop = FALSE]
  yhat_obs <- as.vector(x_obs %*% beta_hat)
  yhat_mis <- as.vector(x_mis %*% beta_dot)

  if (is.null(donors)) {
    donors <- round(length(yhat_obs) / 600 + 7)
  }
  donors <- max(1L, min(as.integer(donors), length(yhat_obs)))
  donor_pos <- mice::matchindex(yhat_obs, yhat_mis, donors)

  sigma_hat <- as.numeric(parm$sigma)
  residual_obs <- y[ry] - as.vector(x_obs %*% drop(parm$coef))
  score_obs <- x_obs * (residual_obs / sigma_hat^2)
  information <- -crossprod(x_obs) / (sigma_hat^2 * n)
  d_obs <- t(-solve(information, t(score_obs)))

  pmm_score <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_d <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  pmm_d[ry, ] <- d_obs
  donor_id <- rep(NA_integer_, n)
  donor_id[wy] <- row_id[ry][donor_pos]

  model$setup <- list(
    method = "pmmrw",
    n = sum(ry),
    task = task,
    donors = donors,
    matchtype = matchtype,
    ridge = ridge,
    experimental_kappa = "direct_joint_integrated_downstream"
  )
  model$beta.hat <- beta_hat
  model$beta.dot <- beta_dot
  model$sigma.dot <- sigma_hat
  model$xnames <- colnames(x)
  model$donor_id <- donor_id
  model$pmm_score <- pmm_score
  model$pmm_d <- pmm_d

  y[ry][donor_pos]
}
