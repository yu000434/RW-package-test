\name{giganti_data}
\alias{giganti_data}
\docType{data}
\title{Giganti and Shepherd simulation data}
\description{
A simulated two-phase dataset based on Giganti and Shepherd (2020),
not participant data from their study. Validated values of \code{A} and
\code{D} are available for 1,000 of the 4,000 observations.
}
\usage{giganti_data}
\format{
A data frame with 4,000 rows and 8 variables:
\describe{
  \item{ID}{Observation identifier.}
  \item{X1, X2}{Fully observed continuous covariates.}
  \item{A.star}{Fully observed error-prone measure of the continuous exposure.}
  \item{D.star}{Fully observed error-prone binary outcome, coded 0/1.}
  \item{A}{Validated continuous exposure; missing outside the validation sample.}
  \item{D}{Validated binary outcome, a factor with levels 0 and 1;
    missing outside the validation sample.}
  \item{Intercept}{A constant column equal to 1.}
}
}
\references{
Giganti MJ, Shepherd BE (2020).
Multiple-Imputation Variance Estimation in Studies With Missing or
Misclassified Inclusion Criteria.
\emph{American Journal of Epidemiology}, 189(12), 1628--1632.
\doi{10.1093/aje/kwaa153}.
}
\keyword{datasets}
