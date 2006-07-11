\name{plot2d}
\alias{plot2d}
\title{Plots two dimensional summary of an ideal object}
\description{
  Two dimensional scatterplots of estimated ideal points provided by an \code{ideal} object
}
\usage{
plot2d(x, d1=1, d2=2, start=rownames(x$x)[1],
       overlayCuttingPlanes=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{\link{ideal}}.}
  \item{d1}{integer, the number of the first dimension to plot.  This dimension
    will appear on the horizontal (x) axis.}
  \item{d2}{integer, the number of the second dimension to plot.  This dimension
    will appear on the vertical (y) axis.}
  \item{start}{the number of the iteration with which to begin when
    calculating the posterior mean and confidence interval.}
  \item{overlayCuttingPlanes}{logical, if \code{TRUE}, overlay the
  estimated bill-specific cutting planes}
  \item{...}{other parameters to be passed through to plotting functions.}
}

\details{\code{\link{plot2d}} plots the posterior mean ideal points for each
  legislator in the first dimension against the ones for the second
  dimension.  If party information is available (i.e., as a variable in
  the \code{legis.data} component of the \code{rollcall} object used to
  create the \code{ideal} object), legislators from different parties
  are plotted in different colors with different symbols.

  The usual warning about identification applies \emph{a fortiori} in the
  case of \code{d} > 1.

  The \code{overlayCuttingPlanes} option will overlay estimates of the
  cutting planes separating the Yea position from the Nay position for
  each roll call.  These planes correspond to the locus of points where
  legislators are estimated to be indifferent between a Yea vote and a
  Nay vote (n.b., for models with \code{d > 1} dimensions, the location
  of the Yea and Nay alternatives are not identified, but the directions
  to those alternatives are identified, as is the hyperplane separating
  the Yea and Nay alternatives.  In two dimensions, the cutting plane
  for rollcall vote \code{j} is a line with intercept .. and slope ...,
  with the estimated mean of the posterior densities for these
  parameters used to actually produce the plotted cutting line. 
}

\seealso{\code{\link{ideal}}, \code{\link{plot1d}}, \code{\link{plot.ideal}}
}

\examples{
data(s109)
cl2 <- constrain.legis(s109,
                       x=list("KENNEDY (D MA)"=c(-1,0),
                         "ENZI (R WY)"=c(1,0),
                         "CHAFEE (R RI)"=c(0,-.5)),
                       d=2)

\dontrun{
## too long for examples
id2Constrained <- ideal(s109,
                        d=2,
                        priors=cl2,      ## priors (w constraints)
                        startvals=cl2,   ## start value (w constraints)
                        store.item=TRUE,
                        maxiter=5000,
                        burnin=500,
                        thin=25)
}
## short run for examples
id2Constrained <- ideal(s109,
                        d=2,
                        priors=cl2,      ## priors (w constraints)
                        startvals=cl2,   ## start value (w constraints)
                        store.item=TRUE,
                        maxiter=500,
                        burnin=100,
                        thin=10)
plot(id2Constrained)
plot(id2Constrained,
     overlayCuttingPlanes=TRUE)
}

\keyword{hplot}

