\name{summary.hurdle}
\alias{summary.hurdle}
\alias{print.summary.hurdle}
\alias{coef.hurdle}

\title{summarize hurdle regression models for count data}

\description{
  Provides \code{summary} methods for hurlde regression models for count
  data.
}

\usage{
\method{summary}{hurdle}(object,...)

\method{print}{summary.hurdle}(x, digits = max(3,getOption("digits")-3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object inheriting from class \code{"hurdle"}}
  \item{x}{an object of class \code{"summary.hurdle"}, usually a result
    of a call to \code{summary.hurdle}}
  \item{digits}{the number of significant digits to use when printing}
  \item{...}{further arguments passed to or from other methods}
}
\details{
 \code{print.summary.hurdle} tries to be smart about formatting the
 display of the MLEs, standard errors, etc, essentially using the same
 code as appears in \code{link{print.summary.lm}}
}
\value{
  The function \code{summary.hurdle} computes and returns a list of
  summary statistics from the zero-inflated regression model, including

  \item{coefficents}{a matrix, with columns for the MLEs, their standard
    errors, \eqn{z}{z}-statistic, and corresponding (two-sided)
    \eqn{p}{p}-value.}
  \item{vc}{The estimated variance-covariance matrix of the MLEs}
  \item{beta}{The MLEs from the count component of the model}
  \item{gamma}{The MLEs from the hurdle component of the model}
  \item{theta}{If a negative binomial count model is fit, the
    MLE of the over-dispersion parameter}
  \item{llh}{The value of the log-likelihood function at the MLEs}
}
\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{hurdle}}}
\examples{data(bioChemists)
hp <- hurdle(count=art ~ .,
              x = ~ fem + mar + kid5 + phd + ment,
              z = ~ fem + mar + kid5 + phd + ment,
              data=bioChemists,trace=TRUE)
summary(hp)

hnb <- hurdle(count=art ~ .,
              x = ~ fem + mar + kid5 + phd + ment,
              z = ~ fem + mar + kid5 + phd + ment,
              dist="negbin",
              data=bioChemists,trace=TRUE)
summary(hnb)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
