\name{vuong}
\alias{vuong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Vuong's non-nested hypothesis test}
\description{
  compares two models fit to the same data that do not nest.
}
\usage{
vuong(m1, m2, digits = getOption("digits"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{model 1, an object inheriting from class \code{glm},
    \code{negbin} or \code{zeroinfl}}
  \item{m2}{model 2, as for model 1}
  \item{digits}{significant digits in printed result}
}
\details{
  The Vuong non-nested test is based on a comparison of the predicted
  probabilities of two models that do not nest.  Examples include
  comparisons of zero-inflated count models with their non-zero-inflated
  analogs (e.g., zero-inflated Poisson versus ordinary Poisson, or
  zero-inflated negative-binomial versus ordinary negative-binomial).  A
  large, positive test statistic provides evidence of
  the superiority of model 1 over model 2, while a large, negative test statistic is
  evidence of the superiority of model 2 over model 1.  Under the null
  that the models are indistinguishable, the test
  statistic is asymptotically distributed standard normal.

  The function will fail if the models do not contain identical values
  in their respective components named \code{y} (the value of the
  response being modeled).
}
\value{
  nothing returned, prints the test-statistic and \eqn{p} value and
  exits silently.
}
\references{Vuong, Q.H. 1989. "Likelihood ratio tests for model
  selection and non-nested hypotheses." \emph{Econometrica}. 57:307-333.}
\author{Simon Jackman \email{jackman@stanford.edu}}
\examples{
data(bioChemists)
glm1 <- glm(art ~ .,
            data=bioChemists,
            family=poisson,trace=TRUE)  ## poisson GLM

zip <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,
                z = ~ fem + mar + kid5 + phd + ment,
                dist="poisson",        
                data=bioChemists,trace=TRUE) ## zero-infl poisson

vuong(glm1,zip)

## compare negbin with zero-inflated negbin
require(MASS)
nb1 <- glm.nb(art ~ .,
              data=bioChemists,
              trace=TRUE)
zinb <- zeroinfl(count=art ~ .,
                 x = ~ fem + mar + kid5 + phd + ment,
                 z = ~ fem + mar + kid5 + phd + ment,
                 dist="negbin",
                 data=bioChemists,trace=TRUE)
vuong(nb1,zinb)
}
\keyword{models}% at least one, from doc/KEYWORDS
