#' Fit an RW analysis to multiply imputed data
#'
#' @param data A `mids` object created with `mice(..., tasks = "train")`.
#' @param expr An `lm` or binomial `glm` expression.
#'
#' @return An `rw_fit` object.
#' @export
with_rw <- function(data, expr) {
  if (!inherits(data, "mids") || is.null(data$models)) {
    stop("`data` must be a mids object created with `tasks = \"train\"`.")
  }
  missing <- lapply(data$data[names(data$models)], is.na)
  variables <- names(missing)[vapply(missing, any, logical(1))]
  methods <- unname(data$method[variables])
  if (any(!methods %in% c("norm", "logreg", "pmmrw"))) {
    stop("Supported imputation methods are norm, logreg, and pmmrw.")
  }
  if (sum(methods == "pmmrw") > 1L) {
    stop("The current PMM implementation supports one pmmrw variable.")
  }

  expression <- substitute(expr)
  environment <- parent.frame()

  results <- lapply(seq_len(data$m), function(p) {
    completed <- mice::complete(data, p)
    rownames(completed) <- seq_len(nrow(completed))

    imputation <- lapply(variables, function(variable) {
      model <- data$models[[variable]][[p]]
      if (data$method[[variable]] == "pmmrw") {
        pmm_component(model)
      } else {
        parametric_component(completed, model, variable, missing[[variable]])
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

rw_variance <- function(object, pmm_kappa = NULL, pmm_columns = NULL, donor_id = NULL) {
  results <- object$results
  m <- object$m
  n <- object$n
  if (is.null(donor_id)) {
    u_sum <- Reduce(`+`, lapply(results, `[[`, "U"))
  } else {
    u_sum <- Reduce(`+`, lapply(seq_len(m), function(p) {
      U <- results[[p]]$U
      for (i in which(!is.na(donor_id[, p]))) {
        donor <- donor_id[i, p]
        U[donor, ] <- U[donor, ] + U[i, ]
        U[i, ] <- 0
      }
      U
    }))
  }
  u_bar <- u_sum / m
  omega <- crossprod(u_bar) / n

  kappa_sum <- alpha_sum <- d_bar_sum <- 0
  for (p in seq_len(m)) {
    kappa_sum <- kappa_sum + crossprod(results[[p]]$U, results[[p]]$S_mis_imp)
    alpha_sum <- alpha_sum + crossprod(results[[p]]$d)
    d_bar_sum <- d_bar_sum + results[[p]]$d
  }
  kappa <- kappa_sum / (n * m)
  if (!is.null(pmm_kappa)) {
    if (!identical(dim(pmm_kappa), c(nrow(kappa), length(pmm_columns)))) {
      stop("`pmm_kappa` has incompatible dimensions.")
    }
    kappa[, pmm_columns] <- pmm_kappa
  }
  alpha <- alpha_sum / (n * m)
  d_bar <- d_bar_sum / m
  correction <- kappa %*% t(d_bar) %*% u_bar
  delta <- omega + kappa %*% alpha %*% t(kappa) + (correction + t(correction)) / n
  tau <- Reduce(`+`, lapply(results, `[[`, "tau")) / (m * n)
  tau_inv <- solve(tau)
  tau_inv %*% delta %*% t(tau_inv) / n
}

#' Pool an RW analysis
#'
#' @param object An `rw_fit` object returned by [with_rw()].
#' @param pmm_kappa Optional PMM cross-term matrix.
#'
#' @return An `rw_pool` object.
#' @export
pool_rw <- function(object, pmm_kappa = NULL) {
  if (!inherits(object, "rw_fit")) stop("`object` must be returned by with_rw().")
  methods <- object$mids$method[names(object$mids$models)]
  pmm_variable <- names(methods)[methods == "pmmrw"]
  pmm_variable <- pmm_variable[vapply(pmm_variable, function(x) {
    anyNA(object$mids$data[[x]])
  }, logical(1))]
  donor_id <- columns <- NULL
  if (length(pmm_variable) == 1L) {
    donor_id <- extract_donor_id(object$mids, pmm_variable)
    columns <- pmm_columns(object$mids, pmm_variable)
    if (is.null(pmm_kappa)) {
      pmm_kappa <- pmm_kappa(object, pmm_variable, pmm_score_lm)
    }
  } else if (!is.null(pmm_kappa)) {
    stop("`pmm_kappa` requires one pmmrw variable.")
  }
  template <- coef(object$results[[1L]]$model)
  estimates <- vapply(object$results, function(x) coef(x$model), numeric(length(template)))
  if (is.null(dim(estimates))) estimates <- matrix(estimates, nrow = length(template))
  estimate <- rowMeans(estimates)
  names(estimate) <- names(template)
  variance <- rw_variance(object, pmm_kappa, columns, donor_id)
  dimnames(variance) <- list(names(estimate), names(estimate))
  structure(list(estimate = estimate, variance = variance, m = object$m,
                 n = object$n, model = object$results[[1L]]$model,
                 call = match.call()), class = "rw_pool")
}

#' @export
coef.rw_pool <- function(object, ...) object$estimate

#' @export
vcov.rw_pool <- function(object, ...) object$variance

#' @export
summary.rw_pool <- function(object, ...) {
  se <- sqrt(diag(object$variance))
  statistic <- object$estimate / se
  if (inherits(object$model, "glm")) {
    critical <- qnorm(0.975)
    p_value <- 2 * pnorm(-abs(statistic))
  } else {
    df <- df.residual(object$model)
    critical <- qt(0.975, df)
    p_value <- 2 * pt(-abs(statistic), df)
  }
  data.frame(term = names(object$estimate), estimate = object$estimate,
             std.error = se, statistic = statistic, p.value = p_value,
             conf.low = object$estimate - critical * se,
             conf.high = object$estimate + critical * se,
             row.names = NULL)
}

#' @export
print.rw_pool <- function(x, ...) {
  cat("Robins-Wang pooled results\n")
  print(summary(x), row.names = FALSE)
  invisible(x)
}

imputation_design <- function(data, model) {
  # Use the same factor coding as MICE, in the recorded coefficient order.
  x <- utils::getFromNamespace("obtain.design", "mice")(
    data, stats::as.formula(model$formula))
  columns <- model$xnames
  columns[columns == ""] <- "(Intercept)"
  x[, columns, drop = FALSE]
}
