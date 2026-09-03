# Evaluates the analysis and imputation components in each completed dataset.

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

    model <- eval(expression, data, environment)
    analysis <- analysis_component(model, nrow(data))
    imputation <- lapply(variables, function(variable) {
      stored_model <- imps$models[[variable]][[p]]
      if (imps$method[[variable]] == "pmmrw") {
        pmm_component(stored_model)
      } else {
        parametric_component(data, stored_model, variable)
      }
    })

    c(list(model = model), analysis,
      list(S_mis_imp = do.call(cbind, lapply(imputation, `[[`, "S_mis_imp")),
           d = do.call(cbind, lapply(imputation, `[[`, "d"))))
  })

  list(results = results, m = imps$m, n = nrow(imps$data), mids = imps)
}
