\name{rw}
\alias{with_rw}
\alias{pool_rw}
\title{Fit and pool an RW analysis}
\description{
Fit an analysis to each completed dataset and combine the estimates using
Robins-Wang variance estimation.
}
\usage{
with_rw(data, expr)
pool_rw(object, pmm_kappa = NULL)
}
\arguments{
\item{data}{A \code{mids} object from \code{mice(..., tasks = "train")}
using \code{norm}, \code{logreg}, or \code{pmmrw}.}
\item{expr}{An unweighted \code{lm} or binomial \code{glm} expression.
Use a numeric response for \code{lm}; for \code{glm}, use the default logit
link and a 0/1 or two-level factor response.}
\item{object}{An \code{rw_fit} object from \code{with_rw()}.}
\item{pmm_kappa}{Optional PMM cross-term matrix.}
}
\value{
\code{with_rw()} returns an \code{rw_fit} with the fitted models, RW components
and original imputations. \code{pool_rw()} returns an \code{rw_pool} with pooled
estimates and variance. Use \code{coef()}, \code{vcov()} and \code{summary()}
to extract results.
}
\details{
\code{with_rw()} retains completed values and variable types. Rows excluded
by the analysis subset receive zero analysis scores.

Parametric cross terms are calculated from the stored analysis and imputation
scores. PMM pooling supports one numeric PMM variable with fully observed
matching-model predictors and accounts for repeated donor use.

Automatic PMM pooling requires the untransformed PMM variable as the
\code{lm()} response. Analysis predictors and any subset must depend only on
fully observed variables. Other analyses need a matrix from
\code{\link{pmm_kappa}} or \code{\link{pmm_kappa_binomial}}.

\code{summary()} reports standard errors, tests and 95\% confidence intervals.
It uses the first fitted model's residual degrees of freedom for \code{lm()}
and the normal distribution for binomial \code{glm()}.
}
