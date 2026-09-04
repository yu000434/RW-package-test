rw_variance <- function(object) {
  results <- object$results
  m <- object$m
  n <- object$n
  u_bar <- Reduce(`+`, lapply(results, `[[`, "U")) / m
  omega <- crossprod(u_bar) / n

  kappa_sum <- alpha_sum <- d_bar_sum <- 0
  for (p in seq_len(m)) {
    kappa_sum <- kappa_sum + crossprod(results[[p]]$U, results[[p]]$S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(results[[p]]$d)
    d_bar_sum <- d_bar_sum + results[[p]]$d
  }
  kappa <- kappa_sum / (n * m)
  alpha <- alpha_sum / (n * m)
  d_bar <- d_bar_sum / m
  correction <- kappa %*% t(d_bar) %*% u_bar
  delta <- omega + kappa %*% alpha %*% t(kappa) + (correction + t(correction)) / n
  tau <- Reduce(`+`, lapply(results, `[[`, "tau")) / (m * n)
  tau_inv <- solve(tau)
  tau_inv %*% delta %*% t(tau_inv) / n
}

#' Pool an RW analysis
#'
#' @param object An `rw_fit` object returned by [with_rw()].
#'
#' @return An `rw_pool` object containing pooled estimates and RW variance.
#' @export
pool_rw <- function(object) {
  if (!inherits(object, "rw_fit")) stop("`object` must be returned by with_rw().")
  template <- coef(object$results[[1L]]$model)
  estimates <- vapply(object$results, function(x) coef(x$model), numeric(length(template)))
  if (is.null(dim(estimates))) estimates <- matrix(estimates, nrow = length(template))
  estimate <- rowMeans(estimates)
  names(estimate) <- names(template)
  variance <- rw_variance(object)
  dimnames(variance) <- list(names(estimate), names(estimate))
  structure(list(estimate = estimate, variance = variance, m = object$m,
                 n = object$n, model = object$results[[1L]]$model,
                 call = match.call()), class = "rw_pool")
}
