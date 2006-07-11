\name{plot.ideal}
\alias{plot.ideal}
\title{plots an ideal object.}
\description{
  plot of the results of an ideal point estimation contained
  in an object of class \code{ideal}
}
\usage{
plot.ideal(x, conf.int=0.95, start=rownames(x$x)[1], ...)
}

\arguments{
  \item{x}{an object of class \code{ideal}.}
  \item{conf.int}{for "ideal" objects with 1 dimension estimated, the
    level of the confidence interval to plot around the
    posterior mean for each legislator.  If 2 or more dimensions were
    estimated, \code{conf.int} is ignored.}
  \item{start}{the number of the iteration with which to begin when
    calculating the posterior mean and confidence interval.}
  \item{...}{other parameters to be passed through to plotting functions.}
}
\details{
  If the \code{ideal} object comes from fitting a \code{d=1} dimensional
  model, then \code{plot.ideal()}
  produces a plot of the posterior mean of the ideal point estimations for each
  legislator with a confidence interval.  This is the same graph that
  would be produced by calling \code{\link{plot1d}}.  If there are more
  than 30 legislators, only 30 points on this graph will be labelled.  If the
  \code{ideal} object has \code{d=2} dimensions, then \code{\link{plot2d}} is
  called, which plots the (estimated) mean of the posterior density of
  each legislator (i.e., the ideal point is a point in 2-dimensional
  Euclidean space, and the posterior density for each ideal point is a
  bivariate density). 

  For unidimensional and two-dimensional models, if party information is
  available in the \code{rollcall}
  object contained in the \code{ideal} object, legislators from different
  parties are plotted in different colors.  If the \code{ideal} object
  has more than 2 dimensions, \code{plot.ideal()} produces a matrix of
  plots of the mean ideal points of each dimension against the posterior
  mean ideal points of the other dimensions. 
}
\seealso{\code{\link{ideal}}, \code{\link{plot1d}}, \code{\link{plot2d}}
}
\examples{
data(s109)
id1 <- ideal(s109,
             d=1,
             meanzero=TRUE,
             store.item=TRUE,
             maxiter=1000,   ## short run for examples
             burnin=100,
             thin=10)  

plot(id1)
}
\keyword{hplot}

