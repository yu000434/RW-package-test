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
