#' Fit an analysis model to multiply imputed data
#'
#' @param data A `mids` object created with `mice(..., tasks = "train")` using
#'   `norm`, `logreg`, or `pmmrw` imputation.
#' @param expr An `lm` or binomial `glm` expression.
#'
#' @return An `rw_fit` object containing the fitted models and RW components.
#' @export
with_rw <- function(data, expr) {
  if (!inherits(data, "mids") || is.null(data$models)) {
    stop("`data` must be a mids object created with `tasks = \"train\"`.")
  }
  variables <- names(data$models)
  variables <- variables[vapply(variables, function(x) anyNA(data$data[[x]]), logical(1))]
  methods <- unname(data$method[variables])
  if (any(!methods %in% c("norm", "logreg", "pmmrw"))) {
    stop("Supported imputation methods are norm, logreg, and pmmrw.")
  }
  if (sum(methods == "pmmrw") > 1L) {
    stop("The current PMM implementation supports one pmmrw variable.")
  }

  expression <- substitute(expr)
  environment <- parent.frame()
  missing <- lapply(variables, function(x) is.na(data$data[[x]]))
  names(missing) <- variables

  results <- lapply(seq_len(data$m), function(p) {
    completed <- mice::complete(data, p)
    rownames(completed) <- seq_len(nrow(completed))
    for (variable in variables) {
      if (is.factor(completed[[variable]])) {
        completed[[variable]] <- as.integer(completed[[variable]]) - 1L
      }
      completed[[paste0(".imputed_", variable)]] <- missing[[variable]]
    }

    imputation <- lapply(variables, function(variable) {
      model <- data$models[[variable]][[p]]
      if (data$method[[variable]] == "pmmrw") {
        pmm_component(model)
      } else {
        parametric_component(completed, model, variable)
      }
    })
    model <- eval(expression, completed, environment)
    analysis <- analysis_component(model, nrow(completed))

    c(list(model = model), analysis,
      list(S_mis_imp = do.call(cbind, lapply(imputation, `[[`, "S_mis_imp")),
           d = do.call(cbind, lapply(imputation, `[[`, "d"))))
  })

  structure(list(results = results, m = data$m, n = nrow(data$data), mids = data,
                 call = match.call()), class = "rw_fit")
}
