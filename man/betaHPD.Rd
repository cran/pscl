\name{betaHPD}
\alias{betaHPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute and optionally plot beta HDRs}
\description{Compute and optionally plot highest density regions for the Beta
  distribution}
\usage{
   betaHPD(alpha,beta,p=.95,plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{scalar, first shape parameter of the Beta density.  Must be
    greater than 1, see details}
  \item{beta}{scalar, second shape parameter of the Beta density.  Must be
    greater than 1, see details}
  \item{p}{scalar, content of HPD, must lie between 0 and 1}
  \item{plot}{logical flag, if \code{TRUE} then plot the density and
    show the HDR}
}
\details{The Beta density arises frequently in Bayesian models of
  binary events, rates, and proportions, which take on values in the
  unit interval.  For instance, the Beta density is a conjugate prior
  for the unknown success probability in binomial trials.  With shape
  parameters \eqn{\alpha > 1} and \eqn{\beta > 1}, the Beta density is
  unimodal.

  In general, a highest density region of \eqn{f(\theta)} with content
  \eqn{p \in [0,1]} is a set of regions \eqn{\mathcal{Q}} of the
  parameter space with the following properties:
  \deqn{\int_\mathcal{Q} f(\theta) d\theta = p}
  and
  \deqn{f(\theta) > f(\theta^*) \, \forall\
    \theta \in \mathcal{Q},
    \theta^* \not\in \mathcal{Q}}
  For a continuous, unimodal density defined with respect to a single
  parameter (like the Beta case considered here), a HDR region is a
  contiguous interval in \eqn{[0,1]}.
  
  This function uses numerical methods are used to solve for the
  boundaries of a HDR region for the Beta density, via repeated calls
  the functions \code{\link{dbeta}}, \code{\link{pbeta}} and
  \code{\link{qbeta}}. The function \code{\link{uniroot}} is used to
  find points \eqn{v} and \eqn{w} such that \deqn{f(v) = f(w)} subject
  to the constraint \deqn{\int_v^w f(\theta; \alpha, \beta) d\theta = p}
}
\value{
  If the numerical optimization is successful an vector of length 2,
  containing \eqn{v} and \eqn{w}, defined above.    If the optimization
  fails for whatever reason, a vector of \code{NAs} is returned.
}
\author{Simon Jackman <jackman@stanford.edu>}
\seealso{\code{\link{pbeta}}, \code{\link{qbeta}},
  \code{\link{dbeta}}, \code{\link{uniroot}}}
\examples{
betaHPD(4,5)
betaHPD(2,120)
}
\keyword{distribution}% at least one, from doc/KEYWORDS
