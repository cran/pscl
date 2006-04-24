\name{summary.zeroinfl}
\alias{summary.zeroinfl}
\alias{coef.zeroinfl}
\alias{print.summary.zeroinfl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing Zero-Inflated Regression Models for Count Data}
\description{
  \code{summary} method for class \code{"zeroinfl"}
}
\usage{
  \method{summary}{zeroinfl}(object, ...)

  \method{print}{summary.zeroinfl}(x,
                digits = max(3, getOption("digits") - 3),
                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"zeroinfl"}, usually a result of a
    call to \code{\link{zeroinfl}}}
  \item{x}{an object of class \code{"summary.zeroinfl"}, usually a
    result of a call to \code{summary.zeroinfl}} 
  \item{digits}{the number of significant digits to use when printing}
  \item{...}{further arguments passed to or from other methods}
}
\details{
 \code{print.summary.zeroinfl} tries to be smart about formatting the
 display of the MLEs, standard errors, etc, essentially using the same
 code as appears in \code{link{print.summary.lm}}
}
\value{
  The function \code{summary.zeroinfl} computes and returns a list of
  summary statistics from the zero-inflated regression model, including

  \item{coefficents}{a matrix, with columns for the MLEs, their standard
    errors, \eqn{z}{z}-statistic, and corresponding (two-sided)
    \eqn{p}{p}-value.}
  \item{vc}{The estimated variance-covariance matrix of the MLEs}
  \item{beta}{The MLEs from the count component of the model}
  \item{gamma}{The MLEs from the zero-inflated component of the model}
  \item{theta}{If a zero-inflated negative binomial model is fit, the
    MLE of the over-dispersion parameter}
  \item{llh}{The value of the log-likelihood function at the MLEs}
}
\author{Simon Jackman <jackman@stanford.edu>}
\seealso{\code{\link{zeroinfl}}, \code{\link{predict.zeroinfl}}

  Methods are supplied for the generic functions \code{link{coef}} and
  \code{link{logLik}}, for objects of class \code{"zeroinfl"}.
}
\examples{
data(bioChemists)
zip <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,
                z = ~ fem + mar + kid5 + phd + ment,
                data=bioChemists,trace=TRUE)
summary(zip)
coef(zip)
logLik(zip)
}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{regression}% __ONLY ONE__ keyword per line
