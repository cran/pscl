\name{dropUnanimous.rollcall}
\alias{dropUnanimous.rollcall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{drop unanimous votes from rollcall objects}
\description{
  drop unanimous votes from rollcall objects.
}
\usage{
\method{dropUnanimous}{rollcall}(obj, lop = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{object of class \code{\link{rollcall}}}
  \item{lop}{numeric, non-negative integer, less than number of
    legislators represented in \code{obj}.  Roll calls with \code{lop}
    or fewer legislators voting in the minority are dropped.}
}
\details{
  This function is a wrapper to \code{\link{dropRollCall}} with the
  \code{lop} component of \code{dropList} set to the value of \code{lop}
  specified in \code{dropUnanimous}.  
}
\value{
  A \code{\link{rollcall}} object.}
}
\author{Simon Jackman \email{jackman@stanford.edu}}
\seealso{
  \code{\link{dropUnanimous.matrix}},
  \code{\link{dropRollCall}},
  \code{\link{rollcall}},
  \code{\link{summary.rollcall}},
  \code{\link{ideal}}
}
\examples{
data(s109)

## the following two lines yield an identical rollcall object
s109.working <- dropUnanimous(s109)
s109.working
s109.working <- dropRollCall(s109,dropList=list(lop=0))
s109.working
}
\keyword{manip}
