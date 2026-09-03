# Assembles the RW variance after replacing the PMM columns of kappa
# and clustering estimating functions by their realized donor source.

compute_rw_variance <- function(
    fit,
    pmm_kappa,
    pmm_columns,
    donor_id) {
  results <- fit$results
  m <- fit$m
  n <- fit$n
  normalize_donor_id <- utils::getFromNamespace("normalize_donor_id", "rw")
  cluster_U_by_donor <- utils::getFromNamespace("cluster_U_by_donor", "rw")
  donor_id <- normalize_donor_id(donor_id, m, n)

  u_sum <- Reduce(`+`, lapply(results, `[[`, "U"))
  # Reused donor values share one source contribution in Omega.
  u_omega_sum <- Reduce(`+`, lapply(seq_len(m), function(p) {
    cluster_U_by_donor(results[[p]]$U, donor_id[, p])
  }))
  tau_sum <- Reduce(`+`, lapply(results, `[[`, "tau"))
  u_bar_omega <- u_omega_sum / m
  omega <- crossprod(u_bar_omega) / n

  kappa_sum <- alpha_sum <- d_bar_sum <- 0
  for (p in seq_len(m)) {
    kappa_sum <- kappa_sum + crossprod(
      results[[p]]$U,
      results[[p]]$S_mis_imp
    )
    alpha_sum <- alpha_sum + crossprod(results[[p]]$d)
    d_bar_sum <- d_bar_sum + results[[p]]$d
  }
  kappa <- kappa_sum / (n * m)
  if (!identical(dim(pmm_kappa), c(nrow(kappa), length(pmm_columns)))) {
    stop("The PMM kappa has incompatible dimensions.")
  }
  # Non-PMM columns remain the original RW terms.
  kappa[, pmm_columns] <- pmm_kappa
  alpha <- alpha_sum / (n * m)
  d_bar <- d_bar_sum / m

  correction <- kappa %*% t(d_bar) %*% u_bar_omega
  delta <- omega +
    kappa %*% alpha %*% t(kappa) +
    (correction + t(correction)) / n
  tau <- tau_sum / (m * n)
  tau_inv <- solve(tau)
  variance <- tau_inv %*% delta %*% t(tau_inv) / n

  list(
    variance = variance,
    kappa = kappa,
    alpha = alpha,
    d_bar = d_bar,
    omega = omega,
    u_bar = u_sum / m,
    u_bar_omega = u_bar_omega
  )
}
