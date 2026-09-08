# Extracts the quantities saved in the simulation results.

var_parts <- function(fit, pooled, target, pmm_kappa = NULL, pmm_columns = NULL,
                      donor_id = NULL) {
  # The final variance comes from pool_rw(); these terms are for reporting only.
  results <- fit$results
  u_bar <- Reduce(`+`, lapply(results, `[[`, "U")) / fit$m
  u_bar_omega <- u_bar
  if (!is.null(donor_id)) {
    u_bar_omega <- Reduce(`+`, lapply(seq_len(fit$m), function(p) {
      U <- results[[p]]$U
      for (i in which(!is.na(donor_id[, p]))) {
        donor <- donor_id[i, p]
        U[donor, ] <- U[donor, ] + U[i, ]
        U[i, ] <- 0
      }
      U
    })) / fit$m
  }
  kappa <- Reduce(`+`, lapply(results, function(x) crossprod(x$U, x$S_mis_imp))) /
    (fit$n * fit$m)
  kappa_norm <- sqrt(sum((if (is.null(pmm_kappa)) kappa else pmm_kappa)^2))
  if (!is.null(pmm_kappa)) kappa[, pmm_columns] <- pmm_kappa
  alpha <- Reduce(`+`, lapply(results, function(x) crossprod(x$d))) / (fit$n * fit$m)
  d_bar <- Reduce(`+`, lapply(results, `[[`, "d")) / fit$m
  tau <- Reduce(`+`, lapply(fit$results, `[[`, "tau")) / (fit$m * fit$n)
  tau_inverse <- solve(tau)
  transform <- function(component) tau_inverse %*% component %*% t(tau_inverse) / fit$n
  correction <- kappa %*% t(d_bar) %*% u_bar_omega
  components <- list(
    omega = crossprod(u_bar_omega) / fit$n,
    kappa = kappa %*% alpha %*% t(kappa),
    cross = (correction + t(correction)) / fit$n
  )
  variance <- vcov(pooled)
  index <- match(target, rownames(variance))
  values <- vapply(components, function(component) transform(component)[index, index], numeric(1))
  total <- variance[index, index]
  stopifnot(isTRUE(all.equal(sum(values), total, tolerance = 1e-10)))
  c(total = total, values, u_bar_norm = sqrt(sum(u_bar^2)),
    u_bar_omega_norm = sqrt(sum(u_bar_omega^2)), kappa_norm = kappa_norm)
}

pmm_checks <- function(imps, variable) {
  errors <- vapply(seq_len(imps$m), function(p) {
    matching <- pmm_inputs(imps, variable, p)$matching
    donors <- imps$models[[variable]][[p]]$setup$donors
    c(max_probability_error = max(abs(colSums(matching$probability) - donors)),
      max_derivative_error = max(vapply(matching$derivative, function(x) {
        max(abs(colSums(x)))
      }, numeric(1))),
      max_intercept_error = max(abs(matching$derivative[["(Intercept)"]])))
  }, numeric(3))
  apply(errors, 1L, max)
}

rubin_parts <- function(fit, target) {
  estimates <- vapply(fit$results, function(result) coef(result$model)[[target]], numeric(1))
  within <- vapply(fit$results, function(result) vcov(result$model)[target, target], numeric(1))
  u_bar <- mean(within)
  b <- stats::var(estimates)
  c(estimate = mean(estimates), u_bar = u_bar, b = b,
    total = u_bar + (1 + 1 / length(estimates)) * b)
}

donor_use <- function(donor_id) {
  donor_id <- unlist(donor_id, use.names = FALSE)
  donor_id <- donor_id[!is.na(donor_id)]
  counts <- table(donor_id)
  c(n_imputed = length(donor_id), n_unique_donors = length(counts),
    max_donor_reuse = max(counts))
}

score_correlation <- function(fit) {
  one <- function(u) {
    u <- u[, apply(u, 2L, stats::sd) > 0, drop = FALSE]
    if (ncol(u) < 2L) return(NA_real_)
    correlation <- stats::cor(u)
    stats::median(abs(correlation[upper.tri(correlation)]))
  }
  values <- vapply(fit$results, function(result) one(result$U), numeric(1))
  if (all(is.na(values))) NA_real_ else stats::median(values, na.rm = TRUE)
}
