run_one_rw <- function(seed, scenario, n, m, k = NA_integer_, imputation = "pmm") {
  if (!imputation %in% c("pmm", "parametric")) stop("Unknown imputation model.")
  gen <- make_rw_data(seed, n = n, type = scenario)
  data <- gen$dat
  predictors <- c("X", "A")
  is_pmm <- imputation == "pmm"

  method <- mice::make.method(data)
  method[] <- ""
  method["Z"] <- if (is_pmm) "pmmrw" else "norm"
  predictor_matrix <- mice::make.predictorMatrix(data)
  predictor_matrix[,] <- 0L
  predictor_matrix["Z", predictors] <- 1L

  if (is_pmm) {
    set.seed(seed)
    imps <- mice::mice(
      data, m = m, method = method, predictorMatrix = predictor_matrix,
      tasks = "train", blots = list(Z = list(donors = k)), print = FALSE
    )
  } else {
    imps <- mice::mice(data, m = m, method = method, predictorMatrix = predictor_matrix,
                       tasks = "train", print = FALSE)
  }
  fit <- fit_rw(
    imps,
    lm(Z ~ X - 1, subset = if (scenario == "robins_1") A == 1 else rep(TRUE, length(A)))
  )
  if (is_pmm) {
    pmm_score <- pmm_kappa_rw(imps, fit, data, "Z", predictors, k, scenario)
    donor_id <- pmm_donors(imps, "Z")
    variance <- compute_rw_variance(fit, pmm_score$kappa, pmm_cols(imps, "Z"), donor_id)
    kappa <- pmm_score$kappa
    diagnostics <- pmm_score$diagnostics
    reuse <- donor_use(donor_id)
  } else {
    variance <- compute_rw_variance(fit)
    kappa <- variance$kappa
    diagnostics <- c(max_membership_probability_error = NA_real_,
                     max_membership_derivative_error = NA_real_,
                     max_intercept_derivative = NA_real_)
    reuse <- c(n_imputed = NA_real_, n_unique_donors = NA_real_, max_donor_reuse = NA_real_)
  }
  rb <- var_parts(fit, variance, "X")
  rr <- rubin_parts(fit, "X")
  s_mis <- mean(vapply(fit$results, function(result) sqrt(sum(result$S_mis_imp^2)), numeric(1)))

  data.frame(
    seed = seed, scenario = scenario, imputation = imputation,
    n = n, validated = NA_integer_, m = m, k = k,
    quadrature_order = NA_integer_, beta0 = gen$beta0,
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
    max_quadrature_tail_bound = NA_real_,
    mean_analysis_n = mean(vapply(fit$results, `[[`, integer(1), "n_analysis")),
    median_abs_score_component_corr = score_correlation(fit)
  )
}
