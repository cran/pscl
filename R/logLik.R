logLik.polr <- function(object, ...){
  if(class(object)!="polr")
    stop("logLik.polr only defined for objects of class polr\n")

  tmp <- object
  if(is.null(tmp$model)){
    stop("no model component found in polr object, try refitting\n")
  }

  y <- tmp$model[,1]
  yNum <- as.numeric(y)
  p <- tmp$fitted.values
  n <- length(y)
  logLike <- rep(NA,n)
  weights <- rep(1,n)
  if(!is.null(object$call$weights))
    weights <- object$model[,"(weights)"]
  for(i in 1:n){
    logLike[i] <- log(p[i,yNum[i]])*weights[i]
  }
  sum(logLike)
}


