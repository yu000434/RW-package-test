\name{pmm}
\alias{pmm_kappa}
\alias{pmm_kappa_binomial}
\title{PMM cross terms}
\description{
Combine expected analysis scores with derivatives of donor-selection
probabilities. Matching-model predictors must be fully observed.
}
\usage{
pmm_kappa(object, variable, expected_score, ...)
pmm_kappa_binomial(object, variable, threshold = -Inf, quadrature_order = 24L)
}
\arguments{
\item{object}{An \code{rw_fit} object from \code{\link{with_rw}}.}
\item{variable}{Name of the numeric PMM-imputed variable.}
\item{expected_score}{Function returning expected analysis-score matrices;
see below.}
\item{...}{Additional arguments passed to \code{expected_score}.}
\item{threshold}{Lower analysis threshold. Match the subset
\code{A > threshold}; the default \code{-Inf} uses all rows.}
\item{quadrature_order}{Number of quadrature points, default 24.}
}
\value{
The PMM columns of the RW cross-term matrix. Pass this matrix to
\code{\link{pool_rw}} as \code{pmm_kappa}.
}
\section{Expected scores}{
\code{pmm_kappa()} calls
\code{expected_score(object, variable, p, missing, donor_mean, ...)}.
Here \code{p} is the imputation index, \code{missing} is a length-\code{n}
logical indicator and \code{donor_mean} contains fitted donor predictions.

Return one matrix per analysis coefficient, in coefficient order. Rows are
observed donors and columns are missing recipients, both in original data
order. Each entry is the expected score for that donor-recipient pair,
accounting for the analysis subset and downstream imputation.

The cross term averages scores times probability derivatives over imputations
and divides by the full sample size. Only the selection probabilities are
differentiated, not the expected scores.
}
\section{Binomial analysis}{
\code{pmm_kappa_binomial()} supports an unweighted
\code{glm(D ~ A, family = binomial())} with the default logit link,
optionally using \code{A > threshold}. \code{A} uses \code{pmmrw} and
\code{D} uses \code{logreg}. The logistic imputation model must include
\code{A} as an untransformed additive predictor; its other predictors must
be fully observed.

The score uses observed responses where available and the recorded logistic
model otherwise. Quadrature integrates over the fitted normal donor
distribution; it does not change the donor pool or imputed values.
}
