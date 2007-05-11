\name{zeroinfl}
\alias{zeroinfl}
\alias{print.zeroinfl}

\encoding{latin1}

\title{Zero-inflated Count Data Regression}
\description{
  Fit zero-inflated regression models for count data via maximum likelihood.
}

\usage{
zeroinfl(formula, data, subset, na.action, weights, offset,
  dist = c("poisson", "negbin", "geometric"),
  link = c("logit", "probit", "cloglog", "cauchit", "log"),
  control = zeroinfl.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model.}
  \item{dist}{character specification of count model family (a log link is 
    always used).}
  \item{link}{character specification of link function in the binary
    zero-inflation model (a binomial family is always used).}
  \item{control}{a list of control arguments specified via
    \code{\link{zeroinfl.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{\dots}{arguments passed to \code{\link{zeroinfl.control}} in the
    default setup.}
}

\details{
  Zero-inflated count models are two-component mixture models
  combining a point mass at zero with a proper count distribution.
  Thus, there are two sources of zeros: zeros may come from 
  both the point mass and from the count component. Usually the count model 
  is a poisson or negative binomial regression (with log link). 
  The geometric distribution is a special case of the negative binomial
  with size parameter equal to 1.
  For modeling the unobserved state (zero vs. count), a binary model is used:
  in the simplest case only with an intercept but potentially containing regressors.
  For this zero-inflation model, a binomial model with different links can be
  used, typically logit or probit.
  
  The \code{formula} mainly describes the count data model, i.e., \code{y ~ x1 + x2}
  specifies a count data regression where all zero counts have the same
  probability of belonging to the zero component. This is equivalent to the
  model \code{y ~ x1 + x2 | 1}, making it more explicit that the zero-inflation
  model only has an intercept. Additionally, further regressors can be added
  to the zero-inflation model so that not all zeros have the same probability
  for belonging to the point mass component or to the count component. A typical
  formula is \code{y ~ x1 + x2 | z1 + z2}. The regressors in the zero and
  the count component can be overlapping (or identical).
  
  All parameters are estimated by maximum likelihood using \code{\link[stats]{optim}},
  with control options set in \code{\link{zeroinfl.control}}.
  Starting values can be supplied, estimated by the EM (expectation maximization)
  algorithm, or by \code{\link[stats]{glm.fit}} (the default). The latter corresponds
  to the first iteration of the EM algorithm and initializes the unobserved state
  as \code{y > 0}, i.e., all zeros are in the perfect component and only the non-zero
  counts in the count component. Standard errors are derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. See
  \code{\link{zeroinfl.control}} for details.
  
  The returned fitted model object is of class \code{"zeroinfl"} and is similar
  to fitted \code{"glm"} objects. For elements such as \code{"coefficients"} or
  \code{"terms"} a list is returned with elements for the zero and count component,
  respectively. For details see below.
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zeroinfl"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{fitted}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.matrix}}. See \code{\link{predict.zeroinfl}} for more details
  on all methods.
}

\value{
  An object of class \code{"zeroinfl"}, i.e., a list with components including
  \item{coefficients}{a list with elements \code{"count"} and \code{"zero"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{a list with the output from the \code{optim} call for
    minimizing the negative log-likelihood,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call,}
  \item{weights}{the case weights used,}
  \item{offset}{the offset vector used (if any),}
  \item{n}{number of observations,}
  \item{df.null}{residual degrees of freedom for the null model (= \code{n - 2}),}
  \item{df.residual}{residual degrees of freedom for fitted model,}
  \item{terms}{a list with elements \code{"count"}, \code{"zero"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{theta}{estimate of the additional \eqn{\theta}{theta} parameter of the
    negative binomial model (if a negative binomial regression is used),}
  \item{SE.logtheta}{standard error for \eqn{\log(\theta)}{log(theta)},}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of all coefficients in the model (derived from the
    Hessian of the \code{optim} output),}
  \item{dist}{character string describing the count distribution used,}
  \item{link}{character string describing the link of the zero-inflation model,}
  \item{linkinv}{the inverse link function corresponding to \code{link},}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"count"} and \code{"zero"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}),}
}

\references{
Cameron, A. Colin and Pravin K. Trevedi. 1998. \emph{Regression Analysis of Count 
Data.} New York: Cambridge University Press.

Cameron, A. Colin and Pravin K. Trivedi. 2005. \emph{Microeconometrics: Methods and Applications}.
Cambridge: Cambridge University Press.

Lambert, Diane. 1992. \dQuote{Zero-Inflated Poisson Regression,
with an Application to Defects in Manufacturing.} \emph{Technometrics}.V34(1):1-14

Zeileis, Achim, Christian Kleiber and Simon Jackman 2007.
\dQuote{Regression Models for Count Data in R.} Report 53,
Department of Statistics and Mathematics, Wirtschaftsuniversität Wien, Research Report Series.
\url{http://epub.wu-wien.ac.at/}.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{zeroinfl.control}}, \code{\link[stats]{glm}},
  \code{\link[stats]{glm.fit}}, \code{\link[MASS]{glm.nb}},
  \code{\link{hurdle}}
}

\examples{
## from Long (1997)
data("bioChemists", package = "pscl")

## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois <- glm(art ~ ., data = bioChemists, family = poisson)
fm_qpois <- glm(art ~ ., data = bioChemists, family = quasipoisson)
fm_nb <- glm.nb(art ~ ., data = bioChemists)

## with simple inflation
## (no regressors for 0 component)
fm_zip <- zeroinfl(art ~ ., data = bioChemists)
fm_zinb <- zeroinfl(art ~ ., data = bioChemists, dist = "negbin", EM = TRUE)

## inflation with regressors (choose starting values by EM)
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- zeroinfl(art ~ . | ., data = bioChemists, EM = TRUE)
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin", EM = TRUE)
}

\keyword{regression}
