\name{zeroinfl}
\alias{zeroinfl}
\alias{print.zeroinfl}
\alias{logLik.zeroinfl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{zero-inflated regression models for count data}
\description{
  fit zero-inflated regression models for count data via maximum likelihood
}
\usage{
zeroinfl(count = y ~ ., x = ~1, z = ~1,
        data = list(),
        link = "logit", dist = "poisson",
        method = "BFGS", trace = FALSE, maxit = 50000,
        na.action = na.omit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{count being modeled, passed as the left-hand side of a \code{\link{formula}}; lowest count must be zero for
    zero-inflated models}
  \item{x}{right-hand side of a \code{\link{formula}} giving covariates
    for the count part of the model}
  \item{z}{right-hand side of a \code{\link{formula}} giving covariates
    for the zero-inflated part of the model}
  \item{data}{a data frame}
  \item{link}{link function for zero-inflated part of the model (choices
  are \code{logit} (default) or \code{probit})}
  \item{dist}{type of count model, \code{"poisson"} (default) or \code{"negbin"}}
  \item{method}{method for maximizing the log-likelihood function, only \code{"BFGS"} and
    \code{"Nelder-Mead"} are supported}
  \item{trace}{logical, if \code{TRUE}, display progress of maximization}
  \item{maxit}{maximum number of iterations in maximization}
  \item{na.action}{method for handling missing data, default is \code{\link{na.omit}}}
}

\details{Zero-inflated count models are a type of two-component mixture
  model, with a component for zero counts, and the other component for
  the positive counts.  Poisson or negative-binomial models are used for
  the count component of the model; logit or probit is typically used to
  model the probability of a zero count.  Covariates are used to predict
  each outcome. \code{\link{optim}} is used to
  find maximum likelihood estimates and to compute a Hessian matrix
  after convergence.}

\value{
  an object of class \code{zeroinfl}, i.e., a list with components including
  \item{stval}{start values used in optimzation}
  \item{par}{Maximum likelihood estimates}
  \item{hessian}{Matrix of second derivatives of the log-likelihood
    function evaluated at the MLEs; computed numerically by \code{\link{optim}}}
  \item{llh}{value of the log-likelihood function at the MLEs}
  \item{y}{vector of counts actually fitted (after any screeing of
    missing data)}
  \item{x}{matrix of covariates used in fitting the count model}
  \item{z}{matrix of covariates used in fitting the zero-inflated component}
}
\references{Lambert, Diane.  1992. "Zero-Inflated Poisson Regression,
With an Application to Defects in Manufacturing."  \emph{Technometrics}.V34(1):1-14

Cameron, A. Colin and Pravin K. Trevedi. 1998. \emph{Regression analysis of count 
data.} New York: Cambridge University Press.

Long, J. Scott. 1997. \emph{Regression Models for Categorical and Limited Dependent 
Variables.} Number 7 in Advanced Quantitative Techniques in the Social 
Sciences. Thousand Oaks, California: Sage. 
}
\author{Simon Jackman \email{jackman@stanford.edu}}
\seealso{\code{\link{summary.zeroinfl}}, \code{\link{predict.zeroinfl}},
  \code{\link{hurdle}}, \code{\link[MASS:glm.nb]{glm.nb}}

  Methods are supplied for the generic functions \code{\link{coef}} and
  \code{\link{logLik}}, for objects of class \code{"zeroinfl"}.
}  
\examples{
data(bioChemists)
zip <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,   ## no z model
                data=bioChemists,trace=TRUE)
zip
summary(zip)
zip2 <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,
                z = ~ fem + mar + kid5 + phd + ment,
                data=bioChemists,trace=TRUE)
zip2
summary(zip2)

zinb <- zeroinfl(count=art ~ .,
                 x = ~ fem + mar + kid5 + phd + ment,
                 z = ~ fem + mar + kid5 + phd + ment,
                 dist="negbin",
                 data=bioChemists,trace=TRUE)
zinb
summary(zinb)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
