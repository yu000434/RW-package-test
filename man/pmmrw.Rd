\name{pmmrw}
\alias{mice.impute.pmmrw}
\alias{extract_donor_id}
\title{PMM imputation and donor recording}
\description{
Use the PMM donor draw implemented in \code{mice} and record donor IDs
and fitted matching models.
}
\usage{
mice.impute.pmmrw(y, ry, x, wy = !ry, task = "impute", model = NULL,
  exclude = NULL, ridge = 1e-05, matchtype = 1L, donors = 5L,
  use.matcher = FALSE, mlocal = 1L, ...)
extract_donor_id(data, variable)
}
\arguments{
\item{y, ry, x, wy}{Data and row indicators supplied by \code{mice}.}
\item{task, model}{Recording arguments supplied by \code{mice} with
\code{tasks = "train"}.}
\item{exclude, matchtype, use.matcher, mlocal}{Keep \code{exclude = NULL},
\code{matchtype = 1}, \code{use.matcher = FALSE} and \code{mlocal = 1}.}
\item{ridge}{Ridge parameter used by \code{mice}.}
\item{donors}{Donor pool size, default 5.}
\item{...}{Additional arguments passed to \code{mice}'s normal draw.}
\item{data}{A \code{mids} object containing a \code{pmmrw} imputation.}
\item{variable}{Name of the PMM-imputed variable.}
}
\value{
\code{mice.impute.pmmrw()} returns values copied from the selected donors.
\code{extract_donor_id()} returns an \code{n} by \code{m} matrix: rows follow
the original data and columns follow imputations. Imputed rows contain the
selected donor's original row number; observed rows contain \code{NA}.
}
\details{
Set \code{method = "pmmrw"} for a numeric variable in
\code{mice(..., tasks = "train")}.

With the supported PMM settings described above, using the same random
seed and donor pool produces the same imputed values as \code{mice.impute.pmm()}.
}
\seealso{\code{\link[mice]{mice.impute.pmm}} for PMM arguments.}
