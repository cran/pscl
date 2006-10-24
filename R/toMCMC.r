## convert ideal object to MCMC object

idealToMCMC <- function(object, burnin=NULL){
  if(class(object)!="ideal")
    stop("idealToMCMC only defined for objects of class ideal")
  
  if(is.null(burnin))
    b <- eval(object$call$burnin)
  keep <- checkBurnIn(object,b)

  return(mcmc(data=object$x[keep,-1],
              start=object$x[keep,1][1],
              thin=eval(object$call$thin),
              end=object$x[nrow(object$x),1])
         )
}
