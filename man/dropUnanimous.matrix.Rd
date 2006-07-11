\name{dropUnanimous.matrix}
\alias{dropUnanimous.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{drop unanimous rollcalls from a roll call matrix}
\description{
  Drop unanimous or lop-sided rollcalls from a roll call matrix.
}
\usage{
dropUnanimous.matrix(obj, lop = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a \code{\link{matrix}}, where rows index 
    legislators, and columns index votes}
  \item{lop}{numeric, non-negative integer, less than number of
    legislators represented in \code{object}.  Roll calls with \code{lop}
    or fewer legislators voting in the minority are dropped.}
}
\details{
  The roll call matrix can only contain \code{c(0,1,NA)}, as is created
  by \code{\link{convertCodes}}. Each column of the roll call matrix is
  inspected, and the minimum of the number of \code{0} (\dQuote{Nay})
  entries and the number of \code{1} (\dQuote{Yea}) entries is noted.
  If this number is less than or equal to \code{lop}, the column is
  deleted from the matrix returned by this function.
}
\value{
  A \code{\link{matrix}}, with \code{n} rows (same as in the input
  matrix) but possibly less than the original number of columns.
}
\author{Simon Jackman \email{jackman@stanford.edu}}
\seealso{\code{\link{dropUnanimous}}, \code{\link{convertCodes}}}
\examples{
data(s109)
mat <- convertCodes(s109)
class(mat)
dim(mat)
dim(dropUnanimous(mat))

## compare with operating on the rollcall object directly
## with dropUnanimous.rollcall
## produces an identical matrix
dim(dropUnanimous(s109)$votes)
}
\keyword{manip}
\keyword{methods}
