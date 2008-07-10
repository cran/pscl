\name{presidentialElections}
\alias{presidentialElections}
\docType{data}
\title{elections for U.S. President, 1932-2004, by state}

\description{
  Data from states and the District of Columbia, Democratic share of the
  presidential vote, 1932-2004.}

\usage{data(presidentialElections)}

\format{
  \itemize{
    \item{state}{name of state, character}
    \item{demVote}{percent of the vote for president won by the
      Democratic candidate}
    \item{year}{numeric, integer}
    \item{south}{logical, \code{TRUE} if state is one of the 11 states
      of the former Confederacy}
  }
}

\note{945 observations, unbalanced panel data.  The imbalance
  arises due to the facts that (a) Hawaii and Alaska contribute data
  from 1960 onwards, (b) the District of Columbia contributes data
  from 1964 onward, and (c) Alabama has missing data for 1948 and
  1964.}

\source{David Leip's Atlas of U.S. Presidential Elections \url{http://uselectionsatlas.org}
}
  
\examples{
data(presidentialElections)
require(lattice)
xyplot(demVote ~ year | state,
       panel=panel.lines,
       ylab="Democratic Vote for President (percent)",
       xlab="Year",
       data=presidentialElections,
       scales=list(y=list(cex=.6),x=list(cex=.35)),
       strip=strip.custom(par.strip.text=list(cex=.6)))
}
\keyword{datasets}
