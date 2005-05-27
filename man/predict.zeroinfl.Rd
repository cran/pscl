\name{predict.zeroinfl}
\alias{predict.zeroinfl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate predictions from zero-inflated regression count models}
\description{
  Generate predicted counts from zero-inflated regression models for
  count data, fit with the \code{\link{zeroinfl}} function, producing
  objects of class \code{zeroinfl}.
}
\usage{
\method{predict}{zeroinfl}(object, newdata, se.fit = FALSE, conf = 0.95,
              MC = 1000, type = c("response", "prob"), na.action = na.pass, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{zeroinfl}}
  \item{newdata}{optionally, a data frame in which to look for variables
    with which to predict. If omitted, the fitted linear predictors are used.}
  \item{se.fit}{logical switch indicating if standard errors on
    predicted values are to be computed}
  \item{conf}{a proportion, width of confidence intervals computed if
    \code{se.fit=TRUE}, defaults to .95}
  \item{MC}{number of Monte Carlo iterates for computing standard errors
  and confidence intervals around predicted values}
  \item{type}{the type of prediction required.  The default is
    \code{"response"}, generating predictions on the scale of the
    observed counts.  Choosing \code{type="prob"} generates a matrix of
  predicted probabilities over the range of counts observed in the data,
as documented in \code{\link{predprob.zeroinfl}}; no standard errors or
confidence intervals are
produced for predicted probabilities. The value of the \code{type}
argument can be abbreviated.}
  \item{na.action}{function determining what should be done with missing
    values in newdata. The default is to predict \code{NA}.}
  \item{...}{further arguments passed to or from other methods.}
}
\details{Monte Carlo methods are used to generate standard errors and
  confidence intervals as follows: \code{MC} samples are drawn from a
  multivariate normal distribution centered on the MLEs of the
  parameters and variance-covariance matrix equal to the MLE's
  variance-covariance matrix.  With sampled set of parameter values,
  predicted values are generated.  The standard deviation of the
  predicted values is returned as \code{se.fit} along with the lower and
  upper quantiles corresponding to confidence level \code{conf}.}
\value{
  a list, with components
  \item{yhat}{a vector of predicted counts}
  \item{mu}{predictions from the count part of the model, on the scale
    of the counts}
  \item{phi}{predicted probabilities from the zero-inflated component of
    the model, the probability of zero count}
  \item{se}{if \code{se.fit=T}, the standard errors of the predicted
    counts}
  \item{lower,uppper}{if \code{se.fit=T}, the lower/upper bounds of confidence
    intervals of content \code{conf}}
  \item{prob}{if \code{type="prob"}, a matrix of predicted
    probabilities, with each row containing predicted probabilities over
    the range of counts observed in the data, i.e., \code{ncol(prob) = length(min(y):max(y))}}
}
\author{Simon Jackman <jackman@stanford.edu>}
\note{Variables are first looked for in newdata and then searched for in the usual way (which will include the environment of the formula used in the fit).}
\seealso{\code{\link{zeroinfl}}, \code{\link{predict.glm}}, \code{\link{predprob.zeroinfl}}}
\examples{
data(bioChemists)
zip <- zeroinfl(count=art ~ .,
                x = ~ fem + mar + kid5 + phd + ment,
                z = ~ fem + mar + kid5 + phd + ment,
                data=bioChemists,trace=TRUE)

newdata <- expand.grid(list(fem="Men",
                            mar="Married",
                            kid5=1,
                            phd=3.103,
                            ment=0:77))

yhat <- predict(zip,newdata=newdata,
                se.fit=TRUE,MC=2500)

\dontrun{
plot(x=newdata$ment,
     y=yhat$yhat,
     xlab="Mentor Articles",
     ylab="Predicted Counts",
     ylim=range(zip$y),
     type="n")
polygon(x=c(newdata$ment,rev(newdata$ment)),
        y=c(yhat$lower,rev(yhat$upper)),
        border=FALSE,
        col=gray(.75))
lines(x=newdata$ment,
      y=yhat$yhat,
      lwd=2)
rug(quantile(bioChemists$ment,c(.05,.50,.95)))
title("Predicted Counts and 95% CI,\nZero-Inflated Poisson Model")
}
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
