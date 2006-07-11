\name{tracex}
\alias{tracex}
\title{Trace Plot of Legislators' Ideal Points.}
\description{
  Produces a trace plot of the samples from the posterior density of
  legislators' ideal points.
}
\usage{
tracex(object, legis=NULL, d=1, conf.int=0.95,
       showAll = FALSE, start=rownames(object$x)[1])
}

\arguments{
  \item{object}{an object of class \code{ideal}.}
  \item{legis}{a vector of either the names of legislators (or
  \code{\link[=pmatch]{partial matches}} of the names as given in the
  \code{\link{dimnames}} of \code{object$x}.}
  \item{d}{numeric, either a scalar or a vector of length two, the
  dimension(s) to be traced.} 
  \item{conf.int}{numeric, the level of the confidence interval on the
  posterior mean to be plotted.}
  \item{showAll}{logical, if \code{TRUE} and \code{length(d)==2},
  display traces for all selected legislators' ideal points on
  the one plot.}
  \item{start}{numeric, integer, the number of the iteration with which
  to begin plotting and calculations.}
}

\details{Produces a trace plot showing the history of the MCMC
  iterations for the ideal point of each of the legislators (partially)
  named in \code{legis}.  For \code{d=1}, each trace plot includes a
  trace over iterations, the cumulative mean, a moving average, the
  MCMC-based estimate of the mean of the posterior, and a confidence
  interval (specified by \code{conf.int}) around the mean of the
  posterior (using the estimated \code{\link[=quantile]{quantiles}}) of
  the respective MCMC iterates).  All of these values are calculated
  beginning with the iteration specified by \code{start}.

  When \code{d} is a vector of length two, a 2-dimensional trace plot is
  displayed, with the \code{d[1]} dimension on the horizontal axis, and
  the \code{d[2]} dimension on the vertical axis.}

\seealso{\code{\link{ideal}}; \code{\link{pmatch}} for matching
  legislators' names.}

\examples{
data(s109)
## short run for demo purposes only
id1 <- ideal(s109,
             meanzero=TRUE,
             maxiter=500,burnin=100,thin=10)
tracex(id1,legis="KENNEDY")

## n.b., no such legislator has Horrendous Goblin
tracex(id1,legis=c("KENNEDY","BOXER","KYL","Horrendous Goblin"))

\dontrun{
id2 <- ideal(s109,
             d=2,
             maxiter=5000,   ## short run, and unidentified!
             burnin=500,
             thin=50)
tracex(id2,d=1,legis=c("KENNEDY","BOXER","KYL","Horrendous Goblin"))
tracex(id2,d=2,legis=c("KENNEDY","BOXER","KYL","Horrendous Goblin"))
tracex(id2,d=1:2,
       legis=c("KENNEDY","BOXER","KYL","Horrendous Goblin"))

## partial matching
tracex(id2,d=1:2,
       legis=c("KENN","BOX","BID","SNO","SPEC","MCCA","KYL",
               "Horrendous Goblin"),
       showAll=TRUE)
}
}

\keyword{hplot}

