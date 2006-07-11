## convert ideal object to MCMC object

idealToMCMC <- function(x, start=rownames(x$x)[1]){
  if(class(x)!="ideal")
    stop("idealToMCMC only defined for objects of class ideal")
  start1 <- checkStart(x,start)
  return(mcmc(data=x$x[start1:nrow(x$x),-1],
              start=as.integer(start),
              thin=x$call$thin,
              end=x$x[nrow(x$x),1])
         )
}
