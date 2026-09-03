load_rw_method <- function(sim_root) {
  source(file.path(sim_root, "R", "method", "pmmrw.R"))
  source(file.path(sim_root, "R", "method", "topk.R"))
  source(file.path(sim_root, "R", "method", "rw_variance.R"))
  source(file.path(sim_root, "R", "method", "metrics.R"))
  source(file.path(sim_root, "R", "method", "rb.R"))
  source(file.path(sim_root, "RW_scenario", "PMM", "R", "generate_rw_data.R"))
}

run_one_rw <- function(seed, scenario, n, m, k, sim_root) {
  load_rw_method(sim_root)
  gen <- make_rw_data(seed, n = n, type = scenario)
  data <- gen$dat
  predictors <- c("X", "A")

  method <- mice::make.method(data)
  method[] <- ""
  method["Z"] <- "pmmrw"
  predictor_matrix <- mice::make.predictorMatrix(data)
  predictor_matrix[,] <- 0L
  predictor_matrix["Z", predictors] <- 1L

  set.seed(seed)
  imps <- mice::mice(
    data, m = m, method = method, predictorMatrix = predictor_matrix,
    tasks = "train", blots = list(Z = list(donors = k)), print = FALSE
  )
  fit <- rw::with_rw(
    imps,
    lm(Z ~ X - 1, subset = if (scenario == "robins_1") A == 1 else rep(TRUE, length(A)))
  )
  rb_score <- rb_kappa_rw(
    imps, fit, data, "Z", predictors, k, scenario
  )
  kappa <- rb_score$kappa
  variance <- compute_rw_variance(
    fit, kappa, pmm_cols(imps, "Z"), rw::extract_donor_id(imps, "Z")
  )
  rb <- var_parts(fit, variance, "X")
  rr <- rubin_parts(fit, "X")
  reuse <- donor_use(rw::extract_donor_id(imps, "Z"))
  s_mis <- mean(vapply(fit$results, function(result) {
    sqrt(sum(result$S_mis_imp^2))
  }, numeric(1)))

  data.frame(
    seed = seed, scenario = scenario, n = n, validated = NA_integer_, m = m, k = k,
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
    max_probability_error = rb_score$diagnostics[["max_membership_probability_error"]],
    max_derivative_error = rb_score$diagnostics[["max_membership_derivative_error"]],
    max_intercept_error = rb_score$diagnostics[["max_intercept_derivative"]],
    max_quadrature_tail_bound = NA_real_
  )
}
