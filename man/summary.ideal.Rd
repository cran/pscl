\name{summary.ideal}
\alias{summary.ideal}
\title{Summarize an Ideal object }

\description{
  Provides a summary of the output from ideal point estimation contained
  in an object of class \code{ideal}.
}

\usage{
\method{summary}{ideal}(object, quantiles = c(.025, .975),
              start = rownames(object$x)[1],
              include.beta=FALSE,...)
}

\arguments{

  \item{object}{an object of class \code{\link{ideal}}.}

  \item{quantiles}{a list of quantiles to report for each legislator's
    ideal point and each item's discrimination parameter (if stored in
    the \code{ideal} object).}
  
  \item{start}{the number of the iteration with which to begin when
    computing summary quantities and means.}
  
  \item{include.beta}{whether or not to calculate summary statistics of
    beta, if beta is available.  If the item parameters were not stored
    in the \code{ideal} object, then \code{include.beta} is ignored.}
  
  \item{...}{further arguments passed to or from other functions}
  
}

\value{
  An item of class \code{summary.ideal} with elements:
  \item{object}{the name of the ideal object as an
    \code{\link[=eval]{unevaluated}} \code{\link{expression}}, produced
    by \code{match.call()$object}}
  \item{xResults}{a \code{\link{list}} of length \code{d} (the dimension
    of the fitted model).  Component \code{i} of the list is a
    \code{\link{matrix}} 
    summarizing the MCMC output for the \code{n} legislators' ideal
    points on the \code{i}-th dimension of the model. The columns of
    this matrix contain the mean of the MCMC draws from the posterior
    density of the legislators ideal points, the standard deviation, and
    the requested \code{quantiles}.} 
  \item{bResults}{a \code{\link{list}} of length \code{d+1}, similar
    to \code{xResults}, but containing summaries of the bill parameters;
    i.e., there are \code{d} discrimination parameters per bill, plus an
    intercept.  If the bill/item parameters were not stored when
    \code{\link{ideal}} was called (\code{store.item=FALSE}), or
    \code{include.beta=FALSE}, then \code{bResults}
    is a list of length zero.}
  \item{bSig}{a \code{link{list}} of length \code{d}, each component a
    vector of length \code{m}, of mode \code{logical}, equal to
    \code{TRUE} if the corresponding discrimination parameter is
    distinguishible from zero; see Details.  If \code{store.item} was
    set to \code{FALSE} when \code{\link{ideal}} was invoked, then
    \code{bSig} is a list of length zero.}   
  \item{party.quant}{if party information is available through the
    \code{rollcall} object that was used to run \code{ideal}, then
    \code{party.quant} gives the posterior mean of the legislators'
    ideal points by party, by dimension.  If no party information is
    available, then \code{party.quant=NULL}.}
}
\details{The tests of whether a discrimination parameters are
  distinguishible from zero first checks to see if the two most extreme
  \code{quantiles} are symmetric (e.g., as are the default value of .025
  and .975).  If so, the corresponding quantiles of the MCMC samples for
  each discrimination parameter are inspected to see if they have
  the same sign. If they do, then the corresponding discrimination
  parameter is flagged as distinguishible from zero; otherwise not.
}
\seealso{\code{\link{ideal}}}
\author{Simon Jackman \email{jackman@stanford.edu}}
\examples{
## fake example
set.seed(314159265)
fakeData <- matrix(sample(x=c(0,1),size=1000,replace=TRUE),
                   10,100)
rc <- rollcall(fakeData)
## short-run for demo purposes
idFake <- ideal(rc,maxiter=1000,burnin=100,thin=10)
summary(idFake)

## Supreme Court Example
data(sc9497)
rc <- rollcall(data=sc9497$votes,
               legis.names=sc9497$legis.names,
               desc=sc9497$desc)
id1 <- ideal(rc)
summary(id1)

\dontrun{
data(s109)
cl2 <- constrain.legis(s109,
                       x=list("KENNEDY (D MA)"=c(-1,0),
                         "ENZI (R WY)"=c(1,0),
                         "CHAFEE (R RI)"=c(0,-.5)),
                       d=2)
id2Constrained <- ideal(s109,
                        d=2,
                        priors=cl2,      ## priors (w constraints)
                        startvals=cl2,   ## start value (w constraints)
                        store.item=TRUE,
                        maxiter=5000,
                        burnin=500,
                        thin=25)

summary(id2Constrained,
        include.items=TRUE)
}
}

\keyword{classes}

