\name{s109}
\alias{s109}
\docType{data}
\title{rollcall object, 109th U.S. Senate}
\description{
  A sample rollcall object, generated using Jeff Lewis' quasi-real-time
  collection of the rollcalls of the 109th U.S. Senate (2005-2006).
}

\usage{data(s109)}

\format{
  A \code{\link{rollcall}} object containing the recorded votes of the 109th
  U.S. Senate, plus information identifying the legislators and the
  rollcalls.
}

\details{
  The most recent rollcall in this version of the data set is dated June
  29, 2006.
}

\source{
  Jeff Lewis' web site:
  \url{http://adric.sscnet.ucla.edu/rollcall/static/S109.ord}
  
  Information identifying the votes is available at
  \url{http://adric.sscnet.ucla.edu/rollcall/static/S109desc.csv}
}

\examples{
data(s109)
require(pscl)
is(s109,"rollcall")
s109
summary(s109)
summary(s109,verbose=TRUE)
\dontrun{
## how s109 was created
s109 <- readKH("http://adric.sscnet.ucla.edu/rollcall/static/S109.ord",
               desc="109th U.S. Senate",
               debug=TRUE)
url <- "http://adric.sscnet.ucla.edu/rollcall/static/S109desc.csv"
s109$vote.data <- data.frame(read.csv(file=url,header=TRUE))
s109$vote.data$date <- as.Date(s109$vote.data$date,
                               format="\%Y-\%m-\%d")
dimnames(s109$votes)[[2]] <- paste(s109$vote.data$session,
                                   s109$vote.data$number,sep="-")
}
}
\keyword{datasets}
