# Constructs the analysis and imputation components used by the RW variance.

imputation_score <- function(data, model, variable) {
  beta <- setNames(as.numeric(model$beta.dot), model$xnames)
  names(beta)[names(beta) == ""] <- "(Intercept)"
  terms <- setdiff(names(beta), "(Intercept)")
  x <- as.matrix(data[terms])
  if ("(Intercept)" %in% names(beta)) x <- cbind(`(Intercept)` = 1, x)
  x <- x[, names(beta), drop = FALSE]
  y <- data[[variable]]
  imputed <- data[[paste0(".imputed_", variable)]]
  observed <- !imputed

  if (model$setup$method == "logreg") {
    mean <- plogis(drop(x %*% beta))
    score <- x * (y - mean)
    information <- -crossprod(x[observed, , drop = FALSE] *
                                sqrt(mean[observed] * (1 - mean[observed]))) / nrow(data)
  } else {
    sigma2 <- model$sigma.dot^2
    residual <- y - drop(x %*% beta)
    score <- cbind(x * (residual / sigma2),
                   sigma2 = 0.5 * (-1 / sigma2 + residual^2 / sigma2^2))
    x_obs <- x[observed, , drop = FALSE]
    residual_obs <- residual[observed]
    p <- ncol(x)
    information <- matrix(0, p + 1L, p + 1L)
    information[seq_len(p), seq_len(p)] <- -crossprod(x_obs) / sigma2
    information[seq_len(p), p + 1L] <- -drop(crossprod(x_obs, residual_obs)) / sigma2^2
    information[p + 1L, seq_len(p)] <- information[seq_len(p), p + 1L]
    information[p + 1L, p + 1L] <-
      sum(1 / (2 * sigma2^2) - residual_obs^2 / sigma2^3)
    information <- information / nrow(data)
  }

  list(S_mis_imp = score * imputed,
       d = t(-solve(information, t(score * observed))))
}

analysis_score <- function(model, n) {
  x <- model.matrix(model)
  y <- model.response(model.frame(model))
  residual <- y - fitted(model)

  if (inherits(model, "glm")) {
    U <- x * residual
    mu <- predict(model, type = "response")
    link <- predict(model, type = "link")
    weight <- family(model)$mu.eta(link)^2 / family(model)$variance(mu)
    tau <- -crossprod(x * sqrt(weight))
  } else {
    sigma2 <- summary(model)$sigma^2
    U <- x * residual / sigma2
    tau <- -crossprod(x) / sigma2
  }

  full_U <- matrix(0, n, ncol(U), dimnames = list(NULL, names(coef(model))))
  full_U[as.integer(rownames(model.frame(model))), ] <- U
  list(U = full_U, tau = tau, n_analysis = nrow(U))
}

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
    analysis <- analysis_score(model, nrow(data))
    components <- lapply(variables, function(variable) {
      imputation_model <- imps$models[[variable]][[p]]
      if (imps$method[[variable]] == "pmmrw") {
        list(S_mis_imp = imputation_model$pmm_score, d = imputation_model$pmm_d)
      } else {
        imputation_score(data, imputation_model, variable)
      }
    })

    c(list(model = model), analysis,
      list(S_mis_imp = do.call(cbind, lapply(components, `[[`, "S_mis_imp")),
           d = do.call(cbind, lapply(components, `[[`, "d"))))
  })

  list(results = results, m = imps$m, n = nrow(imps$data), mids = imps)
}
