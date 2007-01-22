\name{ideal}
\alias{ideal}
\title{analysis of roll call data (IRT models) via Markov
  chain Monte Carlo methods}

\description{
  Analysis of \code{rollcall} data via the spatial voting model;
  analogous to fitting educational testing data via an item-response
  model.  Model fitting via Markov chain Monte Carlo (MCMC).
}

\usage{
ideal(object, codes = object$codes,
      dropList = list(codes = "notInLegis", lop = 0),
      d = 1, maxiter = 10000, thin = 100, burnin = 5000,
      impute = FALSE,
      normalize = FALSE,
      meanzero = normalize,
      priors = NULL, startvals = "eigen",
      store.item = FALSE, file = NULL,
      verbose=FALSE)
}

\arguments{
  \item{object}{an object of class \code{\link{rollcall}}}
  
  \item{codes}{a \code{\link{list}} describing the types of voting
    decisions in the roll call matrix (the \code{votes} component of the
    \code{\link{rollcall}} \code{object}); defaults to 
    \code{object$codes}, the codes in the rollcall object.} 
  
  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}})
    listing voting decisions, legislators and/or votes to be dropped
    from the analysis; see \code{\link{dropRollCall}} for details.}
  
  \item{d}{numeric, (small) positive integer (defaults to 1).}

  \item{maxiter}{numeric, positive integer, multiple of \code{thin}}

  \item{thin}{numeric, positive integer, thinning interval used for
    recording MCMC iterations.}

  \item{burnin}{number of MCMC iterations to run before recording.  The
    iteration numbered \code{burnin} will be recorded.  Must be a
    multiple of \code{thin}.}

  \item{impute}{\code{\link{logical}}, whether to treat missing entries
    of the rollcall matrix as missing at random, sampling from the
    predictive density of the missing entries at each MCMC iteration.}

  \item{normalize}{\code{\link{logical}}, impose identification with the
    constraint that the ideal points have mean of zero and standard
    deviation one.  This option is only recognized for unidimensional
    models, and is sufficient to locally identify the model parameters
    in this case; more restrictions are required for identification when
    \code{d > 1}.  See Details.}

  \item{meanzero}{to be deprecated/ignored; use \code{normalize} instead.}
  
  \item{priors}{a \code{list} of parameters (means and variances)
    specifying normal priors for the legislators' ideal points. The
    default is \code{NULL}, in which case prior values will be generated
    for both legislators' ideal points and roll call parameters (for the
    ideal points, the default prior parameters are mean zero and
    variance one; for the item parameters the defaults are mean zero and
    variance 100).  If not \code{NULL}, \code{priors} must be a
    \code{list} with as many as four named components \code{xp, xpv, bp,
      bpv}, where
    \item{\code{xp}}{a \code{n} by \code{d} matrix
      of prior \emph{means} for the legislators' ideal points;
      or alternatively, a scalar, which will be replicated to fill a \code{n}
      by \code{d} matrix.}
    \item{\code{xpv}}{a \code{n} by \code{d} matrix of prior
      \emph{precisions} (inverse variances);
      or alternatively, a scalar, which will be replicated to fill a \code{n}
      by \code{d} matrix.} 
    \item{\code{bp}}{a \code{m} by \code{d+1} matrix of prior means for the
      item parameters (with the item difficulty parameter coming last);
      or alternatively, a scalar, which will be replicated to fill a \code{m}
      by \code{d+1} matrix.}
    \item{\code{bpv}}{a \code{m} by \code{d+1} matrix of prior
      precisions for the item parameters;
      or alternatively, a scalar, which will be replicated to fill a \code{m}
      by \code{d+1} matrix.}  
    None of the components should contain \code{NA}.  If any
    of the four possible components are not provided, then the
    corresponding component of \code{priors} is assigned using the default
    values described above.}

  \item{startvals}{either a string naming a method for generating start
    values, valid options are \code{"eigen"} (the default) or
    \code{"random"}; or a \code{list} containing start values for
    legislators' ideal points and item parameters.  See Details.}

  \item{store.item}{\code{\link{logical}}, whether item discrimination
    parameters should be stored.  Storing item discrimination parameters
    can consume a large amount of memory.}

  \item{file}{string, file to write MCMC output.  Default is
    \code{NULL}, in which case MCMC output is stored in memory.  Note
    that post-estimation commands like \code{plot} will not work unless
    MCMC output is stored in memory.}

  \item{verbose}{logical, default
    is \code{FALSE}, which generates relatively little output to the R
    console during execution.}
}

\details{The function fits a \code{d}+1 parameter item-response model to
  the roll call data object, so in one dimension the model reduces
  to the two-parameter item-response model popular in educational testing.
  See References.
  
  \strong{Identification}: The model parameters are \strong{not
    identified} without the user supplying some restrictions on the
  model parameters (translations, rotations and re-scalings of the
  ideal points are observationally equivalent, via offsetting
  transformations of the item parameters).  It is the user's
  responsibility to impose these restrictions; the following brief
    discussion provides some guidance.

  For one-dimensional models, a simple route to identification is the
  \code{normalize} option, which guarantees \emph{local} identification
  (identification up to a 180 rotation of the recovered
  dimension). Near-degenerate\dQuote{spike} priors (priors with
  arbitrarily large precisions) or the \code{constrain.legis} option on
  any two legislators' ideal points ensures \emph{global}
  identification.

  Identification in higher dimensions can be obtained by supplying fixed
  values for \code{d+1} legislators' ideal points, provided the supplied
  points span a \code{d}-dimensional space (e.g., three supplied ideal
  points form a triangle in \code{d=2} dimensions), via the
  \code{constrain.legis} option. In this case the function defaults to
  vague normal priors, but at each iteration the sampled ideal points
  are transformed back into the space of identified parameters, applying
  the linear transformation that maps the \code{d+1} fixed ideal points
  from their sampled values to their fixed values.

  Alternatively, one can impose restrictions on the item parameters via
  \code{constrain.items}.

  Another route to identification is via \emph{post-processing}.  That
  is, the user can run \code{ideal} without any identification
  constraints, but then use the function \code{\link{postProcess}} to
  map the MCMC output from the space of unidentified parameters into the
  subspace of identified parameters.  In fact, the when the
  \code{normalize} option is set to \code{TRUE}, the \code{ideal} object
  is simply post-processed with the \code{normalize} option.
  
  \strong{Start values}.  Start values can be supplied by the user, or
  generated by the function itself. 

  The default method, corresponding to \code{startvals="eigen"}, first
  forms a \code{n}-by-\code{n} correlation matrix from the
  double-centered roll call matrix (subtracting row means, and column
  means, adding in the grand mean), and then extracts the first \code{d}
  principal components (eigenvectors), scaling the eigenvectors by the
  square root of their corresponding eigenvector.  Then, to generate
  start values for the rollcall/item parameters, a series of
  \code{\link[=family]{binomial}} \code{\link[=glm]{glms}} are estimated
  (with a probit \code{\link[=make.link]{link}}), one for each
  rollcall/item, \eqn{j = 1, \ldots, m}.  The votes on the \eqn{j}-th
  rollcall/item are binary responses (presumed to be conditionally
  independent given each legislator's latent preference), and the
  (constrained or unconstrained) start values for legislators are used
  as predictors. The estimated coefficients from these probit models are
  stored to serve as start values for the item discrimination and
  difficulty parameters (with the intercepts from the probit GLMs
  multiplied by -1 so as to make those coefficients difficulty parameters).

  The default \code{eigen} method generates extremely good start values
  for low-dimensional models fit to recent U.S. congresses (where high
  rates of party line voting mean low dimensional models fit well).  The
  \code{eigen} method may be computationally expensive or even
  impossible to implement for \code{rollcall} objects with large numbers
  of legislators.

  The \code{random} method generates start values via iid sampling from
  a N(0,1) density, via \code{\link{dnorm}}, and then uses the probit
  method above to get start values for the rollcall/item parameters.
  
  If \code{startvals} is a \code{list}, it must contain the elements
  \code{xstart} and/or \code{bstart}, which should be matrices.
  \code{xstart} must be of dimensions equal to the number of individuals
  (legislators) by \code{d}.  If supplied, \code{startvals$bstart} must
  be of dimensions number of items (votes) by \code{d}+1.  The
  \code{xstart} and \code{bstart} components cannot contain \code{NA}.
  If \code{xstart} is not supplied when \code{startvals} is a list, then
  start values are generated using the default \code{eiegn} method
  described above, and start values for the rollcall/item parameters are
  regenerated using the probit method, ignoring any user-supplied values
  in \code{startvals$bstart}.  That is, user-supplied values in
  \code{startvals$bstart} are only used when accompanied by a valid set
  of start values for the ideal points in \code{startvals$xstart}.  If 

}
  

\value{a \code{\link{list}} of class \code{ideal} with named components

  \item{n}{\code{\link{numeric}}, integer, number of legislators in the
    analysis, after any subseting via processing the \code{dropList}.}

  \item{m}{\code{\link{numeric}}, integer, number of rollcalls in roll
    call matrix, after any subseting via processing the \code{dropList}.}

  \item{d}{\code{\link{numeric}}, integer, number of dimensions
    fitted.}

  \item{x}{a \code{\link{matrix}} containing the MCMC samples for
    the ideal point of each legislator in each dimension for each
    iteration from \code{burnin} to \code{maxiter}, at an interval of
    \code{thin}.  Rows of the \code{x} matrix index iterations; columns
    index legislators.}

   \item{beta}{a \code{\link{matrix}} containing the MCMC samples
    for the item discrimination parameter for each item in each
    dimension, plus an intercept, for each iteration from \code{burnin}
    to \code{maxiter}, at an interval of \code{thin}. Rows of the
    \code{beta} matrix index MCMC iterations; columns index parameters.}

  \item{xbar}{a \code{\link{matrix}} containing the means of the
    MCMC samples for the ideal point of each legislator in each dimension,
    using iterations \code{burnin} to \code{maxiter}, at an interval of
    \code{thin}; i.e., the column means of \code{x}.}

  \item{betabar}{a \code{\link{matrix}} containing the means of
    the MCMC samples for the vote-specific parameters, using iterations
    \code{burnin} to \code{maxiter}, at an interval of \code{thin};
    i.e., the column means of \code{beta}.}

  \item{call}{an object of class \code{\link{call}}, containing
    the arguments passed to \code{ideal} as unevaluated expressions.}
}

\references{
  Albert, James. 1992. Bayesian Estimation of normal ogive item
  response curves using Gibbs sampling. \emph{Journal of Educational
  Statistics}. 17:251-269.

  Clinton, Joshua, Simon Jackman and Douglas Rivers. 2004. The
  Statistical Analysis of Roll Call Data.  \emph{American Political
  Science Review}.  98:335-370.

  Patz, Richard J. and Brian W. Junker. 1999.  A Straightforward
  Approach to Markov Chain Monte Carlo Methods for Item Response
  Models. \emph{Journal of Education and Behavioral
  Statistics}. 24:146-178.
  
  Rivers, Douglas. 2003.  \dQuote{Identification of Multidimensional
  Item-Response Models.} Typescript.  Department of Political Science,
  Stanford University.
}

\author{Simon Jackman \email{jackman@stanford.edu}, with help from Christina
  Maimone and Alex Tahk.}

\seealso{
  \code{\link{rollcall}}, \code{\link{summary.ideal}},
  \code{\link{plot.ideal}}, \code{\link{predict.ideal}}.
  \code{\link{tracex}} for graphical display of MCMC iterative
  history.

  \code{\link{idealToMCMC}} converts the MCMC iterates in an
  \code{ideal} object to a form that can be used by the \code{coda} library.

  \code{\link{constrain.items}} and
  \code{\link{constrain.legis}} for implementing identifying
  restrictions.

  \code{\link{postProcess}} for imposing identifying restrictions
  \emph{ex post}.
  
  \code{\link[MCMCpack:MCMCirt1d]{MCMCirt1d}} and 
  \code{\link[MCMCpack:MCMCirtKd]{MCMCirtKd}} in the \code{MCMCpack}
  package provide similar functionality to \code{ideal}.
  
}

\examples{
data(s109)

## ridiculously short run for examples
id1 <- ideal(s109,
             d=1,
             normalize=TRUE,
             store.item=TRUE,
             maxiter=500,
             burnin=100,
             thin=10,
             verbose=TRUE)  
summary(id1)

\dontrun{
## more realistic long run
idLong <- ideal(s109,
                d=1,
                priors=list(xpv=1e-12,bpv=1e-12),
                normalize=TRUE,
                store.item=TRUE,
                maxiter=260e3,
                burnin=1e4,
                thin=100)  
}
}
\keyword{models}

