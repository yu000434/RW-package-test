run_one_gs <- function(seed, n, validated, m, k = NA_integer_, imputation = "pmm",
                       threshold = 2L, quadrature_order = 24L) {
  if (!imputation %in% c("pmm", "parametric")) stop("Unknown imputation model.")
  gen <- make_gs_data(seed, obs = n, subsample_n = validated, threshold = threshold)
  data <- gen$dat
  data$D <- factor(data$D, levels = c(0, 1))
  predictors <- c("X1", "X2", "A.star", "D.star")
  is_pmm <- imputation == "pmm"

  method <- mice::make.method(data)
  method[] <- ""
  method[c("A", "D")] <- c(if (is_pmm) "pmmrw" else "norm", "logreg")
  predictor_matrix <- mice::make.predictorMatrix(data)
  predictor_matrix[,] <- 0L
  predictor_matrix["A", predictors] <- 1L
  predictor_matrix["D", c(predictors, "A")] <- 1L

  set.seed(seed)
  if (is_pmm) {
    imps <- mice::mice(
      data, m = m, method = method, predictorMatrix = predictor_matrix,
      tasks = "train", blots = list(A = list(donors = k)), print = FALSE
    )
  } else {
    imps <- mice::mice(data, m = m, method = method, predictorMatrix = predictor_matrix,
                       tasks = "train", print = FALSE)
  }
  fit <- rw::with_rw(imps, glm(D ~ A, family = binomial(), subset = A > threshold))
  if (is_pmm) {
    pmm_score <- pmm_kappa_gs(
      imps, fit, data, "A", predictors, k, threshold = threshold,
      quadrature_order = quadrature_order
    )
    donor_id <- rw::extract_donor_id(imps, "A")
    variance <- compute_rw_variance(fit, pmm_score$kappa, pmm_cols(imps, "A"), donor_id)
    kappa <- pmm_score$kappa
    diagnostics <- pmm_score$diagnostics
    reuse <- donor_use(donor_id)
  } else {
    variance <- compute_rw_variance(fit)
    kappa <- variance$kappa
    diagnostics <- c(max_membership_probability_error = NA_real_,
                     max_membership_derivative_error = NA_real_,
                     max_intercept_derivative = NA_real_, max_quadrature_tail_bound = NA_real_)
    reuse <- c(n_imputed = NA_real_, n_unique_donors = NA_real_, max_donor_reuse = NA_real_)
  }
  rb <- var_parts(fit, variance, "A")
  rr <- rubin_parts(fit, "A")
  s_mis <- mean(vapply(fit$results, function(result) sqrt(sum(result$S_mis_imp^2)), numeric(1)))

  data.frame(
    seed = seed, scenario = "GS", imputation = imputation,
    n = n, validated = validated, m = m, k = k,
    quadrature_order = if (is_pmm) quadrature_order else NA_integer_, beta0 = gen$beta0,
    estimate = rr[["estimate"]], rb_se = sqrt(rb[["total"]]),
    rb_total_var = rb[["total"]], omega_var = rb[["omega"]],
    kappa_var = rb[["kappa"]], cross_var = rb[["cross"]],
    u_bar_norm = sqrt(sum(variance$u_bar^2)),
    u_bar_omega_norm = sqrt(sum(variance$u_bar_omega^2)),
    s_mis_imp_norm = s_mis, kappa_norm = sqrt(sum(kappa^2)),
    rr_se = sqrt(rr[["total"]]),
    rr_total_var = rr[["total"]], rr_u_bar = rr[["u_bar"]], rr_b = rr[["b"]],
    n_imputed = reuse[["n_imputed"]], n_unique_donors = reuse[["n_unique_donors"]],
    max_donor_reuse = reuse[["max_donor_reuse"]],
    max_probability_error = diagnostics[["max_membership_probability_error"]],
    max_derivative_error = diagnostics[["max_membership_derivative_error"]],
    max_intercept_error = diagnostics[["max_intercept_derivative"]],
    max_quadrature_tail_bound = diagnostics[["max_quadrature_tail_bound"]],
    mean_analysis_n = mean(vapply(fit$results, `[[`, integer(1), "n_analysis")),
    median_abs_score_component_corr = score_correlation(fit)
  )
}
