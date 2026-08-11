\name{ms.rlos}
\alias{ms.rlos}
\title{Model Selection for r-Largest Order Statistics}
\usage{
ms.rlos(xdat, sig.ed = 0.1, dmin = 0.03, numh = 8)
}
\arguments{
  \item{xdat}{Numeric matrix or data frame. Rows represent independent blocks and columns contain the largest observations in decreasing order.}
  \item{sig.ed}{Significance level used by the entropy-difference procedures.}
  \item{dmin}{Minimum spacing used when filtering fixed-\eqn{h} rK3D candidates.}
  \item{numh}{Number of initial candidate \eqn{h} values.}
}
\description{
Runs model-selection procedures for r-largest order-statistic data.
}
\details{
The candidate set contains rGLO, rGGD, rGEV, rK4D, rLD, rGD, and fixed-\eqn{h} rK3D models. Selection of \eqn{r} uses entropy-difference procedures, and models evaluated at the same \eqn{r} are compared using AIC.
}
\value{
An object of class \code{msrlos} with components:
\itemize{
  \item \code{surv}: Survival result;
  \item \code{rmed}: r-Median result;
  \item \code{alt.surv}: alternating result initialized by Survival;
  \item \code{alt.rmed}: alternating result initialized by r-Median.
}
}
\author{
Shin Yire, Jihong Park, and Jeong-Soo Park. Correspondence: \email{jspark@jnu.ac.kr}.
}
\references{
Bader B, Yan J, Zhang X (2017). Automated selection of r for the r largest order statistics approach with adjustment for sequential testing. \emph{Statistics and Computing}, 27, 1435--1451.

Shin Y, Park JS (2023). Modeling climate extremes using the four-parameter kappa distribution for r-largest order statistics. \emph{Weather and Climate Extremes}, 39, 100533.

Tawn JA (1988). An extreme-value theory model for dependent observations. \emph{Journal of Hydrology}, 101, 227--250.
}
\examples{
\dontrun{
data(pohang)
set.seed(123)
fit <- ms.rlos(pohang)
fit
fit$surv
}
}
