\name{dropRollCall}
\alias{dropRollCall}
\title{drop user-specified elements from a rollcall object}
\description{
  drop user-specified elements of rollcall object, returning a roll
  call object. 
}
\usage{
dropRollCall(object, dropList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{object}{an object of class \code{\link{rollcall}}}

  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}}) with
    some (or all) of the following components:
    
    \item{\code{codes}}{character or numeric, possibly a vector.  If
      character, it should match the names of \code{object$codes},
      indicating the set of entries in \code{object$votes} to be set to
      \code{NA}.  If numeric, then \code{codes} indicates the entries in
      \code{object$votes} that will be set to \code{NA}.}
    
    \item{\code{lop}}{numeric, non-negative integer, less than number of
      legislators represented in \code{object}.  Roll calls with \code{lop}
      or fewer legislators voting in the minority are dropped.}
    
    \item{\code{legisMin}}{numeric, non-negative integer, less than
      number of roll calls represented in \code{object}.  Legislators with
      \code{legisMin} or fewer votes are dropped.}
    
    \item{\code{dropLegis}}{an \code{\link{expression}} that evaluates to
      mode \code{logical}, vector of length equal to the number of
      legislators represented in \code{object}.  The expression is evaluated
      in the \code{legis.data} component of the rollcall \code{object}.
      Legislators for whom the expression evaluates to \code{TRUE} are
      dropped.}

    \item{\code{dropVotes}}{an \code{\link{expression}} that evaluates to
      mode \code{logical}, vector of length equal to the number of
      rollcalls represented in \code{object}.  The expression is evaluated
      in the \code{votes.data} component of the rollcall \code{object}.
      Rollcalls for which the expression evaluates to \code{TRUE} are
      dropped.}
  }
}
\details{It is often desirable to restrict the analysis of roll call
      data in various ways.  For one thing, unanimous votes provide no
      information discriminating among legislators: hence, a useful
      default for summary and analysis is \code{dropList=list(lop=0)}.
      See the examples for other possibilities, limited only by the
      information supplied in \code{legis.data} and \code{votes.data}.
}
\value{
  An object of class \code{\link{rollcall}} with components modified by
  the subsetting indicated in the \code{dropList}.
  }
\author{Simon Jackman \email{jackman@stanford.edu}}

\note{
  With the exception of \code{codes}, each component of
  \code{dropList} generates a vector of mode \code{\link{logical}},
  either with respect to legislators or votes.  These logical
  vectors are then combined element-wise, such that if any one of
  the subseting restrictions is \code{TRUE} for a particular
  legislator or vote, then that legislator or vote is dropped.  Some
  summaries are reported to the console along the way.

  When \code{dropList} uses the \code{dropLegis} or \code{dropVotes}
  components then \code{dropList} should be constructed as an
  \code{\link{alist}} command; this ensures that the
  \code{dropLegis} and \code{dropVotes} components of
  \code{dropList} are objects of mode \code{\link{expression}}, and
  \code{\link{eval}}uated to mode \code{\link{logical}} in the
  \code{legis.data} and \code{vote.data} \code{\link{environment}}s by
  the function, if possible (rather than being evaluated immediately in
  the environment calling \code{dropRollCall} or constructing
  \code{dropList}).  See the examples.  This is not entirely
  satisfactory, and behavior more like the \code{subset} argument in
  function \code{\link{lm}} would be preferable.
}

\seealso{\code{\link{dropUnanimous}}, \code{\link{summary.rollcall}},
  \code{\link{ideal}}, \code{\link{alist}}.} 

\examples{
data(s109)
s109.working <- dropRollCall(s109,
                             dropList=list(lop=0))
s109.working <- dropRollCall(s109,
                             dropList=list(lop=0,
                               code="notInLegis"))
s109.working <- dropRollCall(s109,
                             dropList=list(lop=3,
                               code="notInLegis"))

## note use of alist, since dropLegis is an expression
dropList <- alist(lop=3,
                 dropLegis=party!="D",
                 code="notInLegis")
s109.working <- dropRollCall(s109,dropList=dropList)

}
\keyword{manip}
