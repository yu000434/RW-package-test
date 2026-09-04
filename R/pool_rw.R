rw_variance <- function(object, pmm_kappa = NULL, pmm_columns = NULL, donor_id = NULL) {
  results <- object$results
  m <- object$m
  n <- object$n
  u_sum <- Reduce(`+`, lapply(results, `[[`, "U"))
  u_omega_sum <- u_sum
  if (!is.null(donor_id)) {
    u_omega_sum <- Reduce(`+`, lapply(seq_len(m), function(p) {
      U <- results[[p]]$U
      for (i in which(!is.na(donor_id[, p]))) {
        donor <- donor_id[i, p]
        U[donor, ] <- U[donor, ] + U[i, ]
        U[i, ] <- 0
      }
      U
    }))
  }
  u_bar <- u_omega_sum / m
  omega <- crossprod(u_bar) / n

  kappa_sum <- alpha_sum <- d_bar_sum <- 0
  for (p in seq_len(m)) {
    kappa_sum <- kappa_sum + crossprod(results[[p]]$U, results[[p]]$S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(results[[p]]$d)
    d_bar_sum <- d_bar_sum + results[[p]]$d
  }
  kappa <- kappa_sum / (n * m)
  if (!is.null(pmm_kappa)) {
    if (!identical(dim(pmm_kappa), c(nrow(kappa), length(pmm_columns)))) {
      stop("`pmm_kappa` has incompatible dimensions.")
    }
    kappa[, pmm_columns] <- pmm_kappa
  }
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
#' @param pmm_kappa Optional PMM cross-term matrix for analyses outside the
#'   automatic linear-model path.
#'
#' @return An `rw_pool` object containing pooled estimates and RW variance.
#' @export
pool_rw <- function(object, pmm_kappa = NULL) {
  if (!inherits(object, "rw_fit")) stop("`object` must be returned by with_rw().")
  methods <- object$mids$method[names(object$mids$models)]
  pmm_variable <- names(methods)[methods == "pmmrw"]
  pmm_variable <- pmm_variable[vapply(pmm_variable, function(x) {
    anyNA(object$mids$data[[x]])
  }, logical(1))]
  donor_id <- columns <- NULL
  if (length(pmm_variable) == 1L) {
    donor_id <- extract_donor_id(object$mids, pmm_variable)
    columns <- pmm_columns(object$mids, pmm_variable)
    if (is.null(pmm_kappa)) pmm_kappa <- pmm_lm_kappa(object, pmm_variable)
  } else if (!is.null(pmm_kappa)) {
    stop("`pmm_kappa` requires one pmmrw variable.")
  }
  template <- coef(object$results[[1L]]$model)
  estimates <- vapply(object$results, function(x) coef(x$model), numeric(length(template)))
  if (is.null(dim(estimates))) estimates <- matrix(estimates, nrow = length(template))
  estimate <- rowMeans(estimates)
  names(estimate) <- names(template)
  variance <- rw_variance(object, pmm_kappa, columns, donor_id)
  dimnames(variance) <- list(names(estimate), names(estimate))
  structure(list(estimate = estimate, variance = variance, m = object$m,
                 n = object$n, model = object$results[[1L]]$model,
                 call = match.call()), class = "rw_pool")
}
