\name{predict.ideal}
\alias{predict.ideal}
\alias{print.predict.ideal}
\title{oredicted orobabilities from an ideal object}
\description{compute predicted probabilities from an \code{ideal}
  object.  This predict method uses the posterior mean values of x and beta to
  make predictions.}
\usage{
\method{predict}{ideal}(object,
                        cutoff=.5,
                        start=rownames(object$x)[1],
                        ...)

\method{print}{predict.ideal}(x,digits=2,...)

}
\arguments{
  \item{object}{an object of class \code{ideal} (produced by
    \code{\link{ideal}}) with item parameters (beta) stored.}
  \item{cutoff}{numeric, a value between 0 and 1, the threshold to be used for sepatating predicted
    probabilities into yea and nay votes.}
  \item{start}{of the recorded MCMC samples, the sample that should be
    used as the first iteration included in prediction calculations.}
  \item{x}{object of class \code{predict.ideal}}
  \item{digits}{number of digits in printed object}
  \item{...}{further arguments passed to or from other methods.}  
}

\details{
  Predicted probabilities are computed using the mean of the posterior
  density of
  of \deqn{x} (ideal points, or latent ability) and \deqn{beta} (bill or
  item parameters).  The percentage correctly predicted
  are determined by counting the percentages of votes with predicted
  probabilities of a Yea vote greater than or equal to the \code{cutoff} as the
  threshold.
}

\value{
  An object of class \code{predict.ideal}, containing:
  \item{pred.probs}{the calculated predicted probability for each
    legislator for each vote.}
  \item{prediction}{the calculated prediction (0 or 1) for each
    legislator for each vote.}
  \item{correct}{for each legislator for each vote, whether the
    prediction was correct.}
  \item{legis.percent}{for each legislator, the percent of votes
    correctly predicted.}
  \item{vote.percent}{for each vote, the percent correctly predicted.}
  \item{yea.percent}{the percent of yea votes correctly predicted.}
  \item{nay.percent}{the percent of nay votes correctly predicted.}
  \item{party.percent}{the average value of the percent correctly
    predicted by legislator, separated by party, if party information
    exists in the \code{rollcall} object used for \code{ideal}.  If no
    party information is available, \code{party.percent = NULL}.}
  \item{overall.percent}{the total percent of votes correctly
    predicted.}
  \item{ideal}{the name of the \code{\link{ideal}} object, which can be
    later \code{\link{eval}}uated}
  \item{desc}{string, the descriptive text from the
  \code{\link{rollcall}} object passed to \code{\link{ideal}}}
}

\seealso{\code{\link{ideal}}, \code{\link{summary.ideal}}, \code{\link{plot.predict.ideal}}}

\examples{
data(s109)

\dontrun{
id1 <- ideal(s109, meanzero=TRUE,
             store.item=TRUE)      ## too long for examples
}

id1 <- ideal(s109,
             d=1,
             meanzero=TRUE,
             store.item=TRUE,  ## need this to be TRUE for predict
             maxiter=1000,
             burnin=100,
             thin=10)  

phat <- predict(id1)
phat         ## print method
}

\keyword{classes}
