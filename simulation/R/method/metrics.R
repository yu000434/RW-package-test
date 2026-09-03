var_parts <- function(fit, variance, target) {
  tau <- Reduce(`+`, lapply(fit$results, `[[`, "tau")) / (fit$m * fit$n)
  tau_inverse <- solve(tau)
  transform <- function(component) {
    tau_inverse %*% component %*% t(tau_inverse) / fit$n
  }
  correction <- variance$kappa %*% t(variance$d_bar) %*% variance$u_bar_omega
  components <- list(
    omega = variance$omega,
    kappa = variance$kappa %*% variance$alpha %*% t(variance$kappa),
    cross = (correction + t(correction)) / fit$n
  )
  index <- match(target, rownames(variance$variance))
  values <- vapply(components, function(component) {
    transform(component)[index, index]
  }, numeric(1))
  c(total = variance$variance[index, index], values)
}

rubin_parts <- function(fit, target) {
  estimates <- vapply(fit$results, function(result) {
    coef(result$model)[[target]]
  }, numeric(1))
  within <- vapply(fit$results, function(result) {
    vcov(result$model)[target, target]
  }, numeric(1))
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
