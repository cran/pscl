\name{predprob.zeroinfl}
\alias{predprob.zeroinfl}

\title{predicted probabilities from zero-inflated regression models}

\description{
  Compute predicted probabilities of counts from a zero-inflated
  regression models.
}

\usage{
\method{predprob}{zeroinfl}(obj,...)
}

\arguments{
  \item{obj}{a fitted object of class inheriting from \code{"zeroinfl"}}
  \item{...}{additional arguments to be passed to \code{\link{predict.zeroinfl}}
}} 

\details{
  Calls \code{\link{predict.zeroinfl}} with the option \code{type="prob"}
}

\value{
  See \code{\link{predict.zeroinfl}}
}

\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{predprob}}, \code{\link{predict.zeroinfl}}}

\examples{
data(bioChemists)
zip <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,
                z = ~ fem + mar + kid5 + phd + ment,
                dist="poisson",        
                data=bioChemists,trace=TRUE) ## zero-infl poisson
phat <- predprob(zip)
newdata <- expand.grid(list(fem="Men",mar="Married",
                            kid5=1,phd=3.103,
                            ment=0:77))
phat <- predprob(zip,newdata=newdata)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
