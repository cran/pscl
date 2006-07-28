\name{hurdle}
\alias{hurdle}
\alias{print.hurdle}
\alias{logLik.hurdle}

\title{hurdle models for count data}

\description{
  Hurdle regression models for count data model, estimation via maximum
  likelihood. 
}

\usage{
hurdle(count = y ~ ., x = ~1, z = ~1, data = list(), link = "logit",
        dist = "poisson", method = "BFGS", trace = FALSE,
        maxit = 50000, na.action = na.omit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{count being modeled, passed as the left-hand side of a \code{link{formula}}; lowest count must be zero}
  \item{x}{right-hand side of a formula giving covariates
    for the count part of the model}
  \item{z}{right-hand side of a formula giving covariates
    for the hurdle part of the model}
  \item{data}{a data frame}
  \item{link}{link function for zero-inflated part of the model (choices
  are \code{logit} or \code{probit})}
  \item{dist}{count model, \code{"poisson"} (default) or \code{"negbin"}}
  \item{method}{method for maximizing the log-likelihood function, only \code{"BFGS"} and
    \code{"Nelder-Mead"} are supported}
  \item{trace}{logical, if \code{TRUE}, display progress of maximization}
  \item{maxit}{maximum number of iterations in maximization}
  \item{na.action}{method for handling missing data, default is \code{\link{na.omit}}}
}
\details{}
\value{An object of class \code{"hurdle"}}
}
\references{}

\author{Simon Jackman \email{jackman@stanford.edu}}

\note{Only minimal support for this class thus far; e.g., no predict method} 

\seealso{\code{\link{zeroinfl}}, \code{\link{summary.hurdle}}}
\examples{
data(bioChemists)
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

\keyword{models}
\keyword{regression}
