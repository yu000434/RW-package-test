run_one_gs <- function(seed, n, validated, m, k, threshold = 2L,
                       quadrature_order = 24L) {
  gen <- make_gs_data(seed, obs = n, subsample_n = validated, threshold = threshold)
  data <- gen$dat
  data$D <- factor(data$D, levels = c(0, 1))
  predictors <- c("X1", "X2", "A.star", "D.star")

  method <- mice::make.method(data)
  method[] <- ""
  method[c("A", "D")] <- c("pmmrw", "logreg")
  predictor_matrix <- mice::make.predictorMatrix(data)
  predictor_matrix[,] <- 0L
  predictor_matrix["A", predictors] <- 1L
  predictor_matrix["D", c(predictors, "A")] <- 1L

  set.seed(seed)
  imps <- mice::mice(
    data, m = m, method = method, predictorMatrix = predictor_matrix,
    tasks = "train", blots = list(A = list(donors = k)), print = FALSE
  )
  fit <- rw::with_rw(imps, glm(D ~ A, family = binomial(), subset = A > threshold))
  pmm_score <- pmm_kappa_gs(
    imps, fit, data, "A", predictors, k, threshold = threshold,
    quadrature_order = quadrature_order
  )
  kappa <- pmm_score$kappa
  donor_id <- rw::extract_donor_id(imps, "A")
  variance <- compute_rw_variance(
    fit, kappa, pmm_cols(imps, "A"), donor_id
  )
  rb <- var_parts(fit, variance, "A")
  rr <- rubin_parts(fit, "A")
  reuse <- donor_use(donor_id)
  s_mis <- mean(vapply(fit$results, function(result) {
    sqrt(sum(result$S_mis_imp^2))
  }, numeric(1)))

  data.frame(
    seed = seed, scenario = "GS", n = n, validated = validated, m = m, k = k,
    quadrature_order = quadrature_order, beta0 = gen$beta0,
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
    max_probability_error = pmm_score$diagnostics[["max_membership_probability_error"]],
    max_derivative_error = pmm_score$diagnostics[["max_membership_derivative_error"]],
    max_intercept_error = pmm_score$diagnostics[["max_intercept_derivative"]],
    max_quadrature_tail_bound = pmm_score$diagnostics[["max_quadrature_tail_bound"]]
  )
}
