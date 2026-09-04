# Collects the model components and assembles the RW variance.

fit_rw <- function(imps, expr) {
  expression <- substitute(expr)
  environment <- parent.frame()
  variables <- names(imps$models)
  variables <- variables[vapply(variables, function(x) anyNA(imps$data[[x]]), logical(1))]
  missing <- lapply(variables, function(x) is.na(imps$data[[x]]))
  names(missing) <- variables

  results <- lapply(seq_len(imps$m), function(p) {
    data <- mice::complete(imps, p)
    for (variable in variables) {
      if (is.factor(data[[variable]])) {
        data[[variable]] <- as.numeric(as.character(data[[variable]]))
      }
      data[[paste0(".imputed_", variable)]] <- missing[[variable]]
    }

    imputation <- lapply(variables, function(variable) {
      stored_model <- imps$models[[variable]][[p]]
      if (imps$method[[variable]] == "pmmrw") {
        pmm_component(stored_model)
      } else {
        parametric_component(data, stored_model, variable)
      }
    })
    model <- eval(expression, data, environment)
    analysis <- analysis_component(model, nrow(data))

    c(list(model = model), analysis,
      list(S_mis_imp = do.call(cbind, lapply(imputation, `[[`, "S_mis_imp")),
           d = do.call(cbind, lapply(imputation, `[[`, "d"))))
  })

  list(results = results, m = imps$m, n = nrow(imps$data), mids = imps)
}

compute_rw_variance <- function(fit, pmm_kappa = NULL, pmm_columns = NULL, donor_id = NULL) {
  results <- fit$results
  m <- fit$m
  n <- fit$n

  u_sum <- Reduce(`+`, lapply(results, `[[`, "U"))
  u_omega_sum <- u_sum
  if (!is.null(donor_id)) {
    # Reused donor values share one source contribution in Omega.
    u_omega_sum <- Reduce(`+`, lapply(seq_len(m), function(p) {
      U <- results[[p]]$U
      for (i in which(!is.na(donor_id[, p]))) {
        j <- donor_id[i, p]
        U[j, ] <- U[j, ] + U[i, ]
        U[i, ] <- 0
      }
      U
    }))
  }
  tau_sum <- Reduce(`+`, lapply(results, `[[`, "tau"))
  u_bar_omega <- u_omega_sum / m
  omega <- crossprod(u_bar_omega) / n

  kappa_sum <- alpha_sum <- d_bar_sum <- 0
  for (p in seq_len(m)) {
    kappa_sum <- kappa_sum + crossprod(results[[p]]$U, results[[p]]$S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(results[[p]]$d)
    d_bar_sum <- d_bar_sum + results[[p]]$d
  }
  kappa <- kappa_sum / (n * m)
  if (!is.null(pmm_kappa)) {
    if (!identical(dim(pmm_kappa), c(nrow(kappa), length(pmm_columns)))) {
      stop("The PMM kappa has incompatible dimensions.")
    }
    kappa[, pmm_columns] <- pmm_kappa
  }
  alpha <- alpha_sum / (n * m)
  d_bar <- d_bar_sum / m

  correction <- kappa %*% t(d_bar) %*% u_bar_omega
  delta <- omega +
    kappa %*% alpha %*% t(kappa) +
    (correction + t(correction)) / n
  tau <- tau_sum / (m * n)
  tau_inv <- solve(tau)
  variance <- tau_inv %*% delta %*% t(tau_inv) / n

  list(variance = variance, kappa = kappa, alpha = alpha, d_bar = d_bar,
       omega = omega, u_bar = u_sum / m, u_bar_omega = u_bar_omega)
}
