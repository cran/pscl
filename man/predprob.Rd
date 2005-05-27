\name{predprob}
\alias{predprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute predicted probabilities from fitted models}
\description{
  Compute predicted probabilities from fitted models, optionally at new
  covariate values
}
\usage{
predprob(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{fitted model object}
  \item{...}{other arguments}
}
\details{
  See documentation for specific methods.
}
\value{
  A matrix of predicted probabilities, each row a vector of
  predicted probabilities over the range of responses seen in the data
  (i.e., \code{min(y):max(y)}), conditional on the values of covariates. 
}
\author{Simon Jackman <jackman@stanford.edu>}
\seealso{\code{\link{predprob.glm}}, \code{\link{predprob.zeroinfl}}}
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
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{models}% at least one, from doc/KEYWORDS
