\name{dropUnanimous}
\alias{dropUnanimous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{drop unanimous votes from rollcall objects and matrices}
\description{
  Drop unanimous votes from rollcall objects and rollcall matrices.
}
\usage{
dropUnanimous(obj, lop = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{object, either of class \code{\link{rollcall}} or
    \code{\link{matrix}}} 
  \item{lop}{numeric, non-negative integer, less than number of
    legislators represented in \code{obj}.  Roll calls with \code{lop}
    or fewer legislators voting in the minority are dropped.  Default is
  0, meaning that unanimous votes are dropped.}
}
\details{
  See documentation for specific methods.  Unanimous votes are the
  equivalent of test items that all subjects score \dQuote{correct} (or
  all subjects scores \dQuote{incorrect}); since there is no variation
  among the legislators/subjects, these votes/items provide no
  information as to latent traits (ideology, preferences, ability).  A
  reasonably large number of rollcalls in any contemporary U.S. Congress
  are unanimous.
}
\value{
  A \code{\link{rollcall}} object or a \code{link{matrix}} depending on
  the class of \code{object}.
}
\author{Simon Jackman \email{jackman@stanford.edu}}
\seealso{\code{\link{dropUnanimous.rollcall}},
  \code{\link{dropUnanimous.matrix}},
  \code{\link{dropRollCall}},
  \code{\link{rollcall}}, \code{\link{summary.rollcall}}, \code{\link{ideal}}}
\examples{
data(s109)
s109.working <- dropUnanimous(s109)
summary(s109.working)
}
\keyword{manip}
