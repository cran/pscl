\name{idealToMCMC}
\alias{idealToMCMC}
\title{Convert and object of class ideal to an MCMC object}
\description{
  Converts the x element of an \code{ideal} object to an MCMC object,
  as used in the \code{coda} package.
}

\usage{
idealToMCMC(x, start=rownames(x$x)[1])
}

\arguments{

  \item{x}{an object of class \code{\link{ideal}}.}

  \item{start}{numeric, integer, the number of the iteration with which
    to begin when converting to an MCMC object.}  }

\value{A \code{\link[coda:mcmc]{mcmc}} object as used by the
  \code{coda} package, starting at
  iteration \code{start}, drawn from the \code{x} component of the
  \code{\link{ideal}} object.}

\seealso{\code{\link{ideal}}, \code{\link[coda:mcmc]{mcmc}}}

\examples{
data(s109)
id1 <- ideal(s109,
             d=1,
             meanzero=TRUE,
             maxiter=1000,   ## short run for demo purposes
             burnin=100,
             thin=10)
id1coda <- idealToMCMC(id1)
summary(id1coda)  
}

\keyword{classes}

