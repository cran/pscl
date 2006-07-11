\name{plot1d}
\alias{plot1d}
\title{Plots 1 dimension of an ideal object}
\description{
  Plots the posterior means of legislator ideal points and
  confidence intervals for a specified dimension of an
  \code{ideal} object.  
}
\usage{
plot1d(x, d=1, conf.int=0.95, start=rownames(x$x)[1],
       showAllNames = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{\link{ideal}}}
  \item{d}{the number of the dimension to plot; only meaningful when
    \code{x} was fit with \code{d} > 1}
  \item{conf.int}{for \code{link{ideal}} objects fit with \code{d=1}
    dimension, the size of the confidence interval to plot around the
    estimated posterior mean for each legislator's ideal point}
  \item{start}{the number of the iteration with which to begin when
    calculating the posterior mean and confidence interval}
  \item{showAllNames}{\code{\link{logical}}, if \code{TRUE}, the
    vertical axis will the names of all legislators. Default is
    \code{FALSE} to reduce clutter on typical-sized graph.}
  \item{...}{other parameters to be passed through to plotting functions}
}
\details{
  \code{plot.ideal()}
  produces a plot of the posterior mean of the ideal point estimations for each
  legislator with a confidence interval.  If there are more
  than 30 legislators, only 30 points on this graph will be labelled.
  If party information is
  available in the \code{rollcall}
  object contained in the \code{ideal} object, legislators from different
  parties are plotted in different colors.  

}

\seealso{\code{\link{ideal}}, \code{\link{plot.ideal}}, \code{\link{plot2d}}}

\examples{
data(s109)
id1 <- ideal(s109,
             d=1,
             meanzero=TRUE,
             store.item=TRUE,
             maxiter=1000,
             burnin=100,
             thin=10)  

plot(id1)  ## same as plot1d(id1) for d=1
plot(id1,showAllNames=TRUE)

\dontrun{
## long run, too long for examples
idLong <- ideal(s109,
                d=1,
                meanzero=TRUE,
                store.item=TRUE,
                maxiter=251e3,
                burnin=1000,
                thin=1e3)  

pdf(file="s109.pdf",
    paper="special",
    height=16,width=8.5)   ## big graphing area
plot(idLong,showAllNames=TRUE)
dev.off()
}
}

\keyword{hplot}

