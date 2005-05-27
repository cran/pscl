\name{odTest}
\alias{odTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test for over-dispersion in count data}
\description{
  Compares the log-likelihoods of a negative binomial regression model
  and a Poisson regression model
}
\usage{
odTest(glmobj, digits = max(3, getOption("digits") - 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{glmobj}{ an object of class negbin produced by glm.nb}
  \item{digits}{ number of digits in printed output}
}

\details{The negative binomial model relaxes the assumption in the
  Poisson model that the (conditional) variance equals the (conditional)
  mean, by estimating one extra parameter.  A likelihood ratio (LR) test
  can be used to test the null hypothesis that the restriction implicit
  in the Poisson model is true. The LR test-statistic has a non-standard
  distribution, even asymptotically, since the negative binomial
  over-dispersion parameter (called theta in glm.nb) is restricted to be
  positive.  The asymptotic distribution of the LR test-statistic has
  probability mass of one half at zero, and a half
  \eqn{\chi^2_1}{chi-square (1)} distribution above zero.  This means
  that if testing at the \eqn{p} = .05 level, we should not reject the
  null unless the LR test-statistic exceeds the critical value
  associated with the \eqn{p} = .025 level; this LR test involves just
  one parameter restriction, so the critical value of the test statistic
  at the \eqn{p} = .05 test statistic is 5.02, instead of the usual 3.8
  (i.e., the .975 quantile of the \eqn{\chi^2_1}{chi-square (1)}
  distribution, versus the .95 quantile).

  A Poisson model is run using \code{\link{glm}} with family set to \code{link{poisson}}, using the
  \code{\link{formula}} in the negbin model object passed as input.  The
  \code{\link{logLik}} functions are used to extract the log-likelihood
  for each model.}

\value{None; prints results and returns silently}
\references{A. Colin Cameron and Pravin K. Trivedi (1998) \emph{Regression
    analysis of count data}. New York: Cambridge University Press.}
\author{Simon Jackman <jackman@stanford.edu>}

\seealso{\code{\link{glm.nb}}, \code{\link{logLik}}}

\examples{
data(bioChemists)
require(MASS)
modelnb <- glm.nb(art ~ .,
                 data=bioChemists,
                 trace=TRUE)
odTest(modelnb)
}
\keyword{regression}% at least one, from doc/KEYWORDS
