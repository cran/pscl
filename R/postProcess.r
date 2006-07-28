## post-process an ideal object
postProcess <- function(object,
                        constraints,
                        debug=FALSE)
{
  if(class(object)!="ideal")
    stop("postProcess only defined for objects of class ideal")

  d <- object$d
  n <- object$n
  m <- object$m
  burnin <- eval(object$call$burnin)
  theIters <- object$x[,1]
  keep <- theIters > burnin
  niter <- dim(object$x)[1]
  
  ## process constraints
  if(!is.list(constraints))
    stop("constraints must be a list")
  nCon <- length(constraints)
  if(nCon != (d+1)){
    cat("postProcess is currently only implements as many constraints\n")
    cat("as there are dimensions plus one.\n")
    stop()
  }
  lengthCon <- lapply(constraints,length)
  if(any(lengthCon != d))
    stop("each constraint must have the same number of dimensions as the fitted model")

  ## form target matrix
  target <- matrix(NA,d+1,d)
  for(i in 1:(d+1))
    target[i,] <- constraints[[i]]

  ## form id vector, where are the named legislators in the ideal object?
  legis.names <- dimnames(eval(object$call$object)$votes)[[1]]
  if(is.null(legis.names)){
    cat("can not find legislator names to match against\n")
    cat(paste("either the original roll call object",
              object$call$object,
              "has been deleted\n"))
    cat("or the vote component of the roll call object has been deleted?\n")
    cat("terminating postProcess with an error\n")
    stop()
  }
 
  legis <- names(constraints)
  ind <- rep(NA,d+1) 
  for(i in 1:nCon){
    p <- grep(pattern=paste("^",legis[i],sep=""),
              x=legis.names)
    if(length(p)==0)
      stop("could not find the named legislator in the rollcall object")
    else
      ind[i] <- p
  }
  cat(paste("matching legislators",legis,"\n"))

  ## initialize output objects
  newX <- NA * object$x
  newX[,1] <- theIters
  dimnames(newX) <- dimnames(object$x)
  
  haveBeta <- eval(object$call$store.item)
  if(haveBeta){
    cat("will also transform item/bill parameters\n")
    newBeta <- NA * object$beta
    newBeta[,1] <- object$beta[,1]
    dimnames(newBeta) <- dimnames(object$beta)
  }

  ## now loop over iterations
  for(iter in 1:niter){
    thisIter <- theIters[iter]
    cat(paste("post-processing iteration",thisIter,"\n"))
    x0 <- matrix(object$x[iter,-1],
                 ncol=d,
                 byrow=TRUE)
    trans <- affineTrans(x=matrix(x0[ind,],nrow=d+1,ncol=d,byrow=FALSE),
                         target=target)
    tmpX <- cbind(x0,1)%*%trans
    newX[iter,-1] <- as.vector(t(tmpX[,1:d]))

    ## now transform beta (and alpha), if available
    if(haveBeta){
      itMat <- try(solve(trans))
      if(!inherits(itMat,"try-error")){
        beta0 <- matrix(object$beta[iter,-1],
                        nrow=d+1,
                        byrow=FALSE)
        tmpBeta <- itMat%*%beta0
        newBeta[iter,-1] <- as.vector(tmpBeta)
      }

      if(debug){
        muPP <- tmpX%*%tmpBeta
        mu <- cbind(x0,1)%*%beta0
        cat("sanity check, comparison of predictions from original and post-processed:\n")
        print(summary(as.vector(mu-muPP)))
      }
    }

  }

  ## clean up objects for return to user as a proper ideal object
  newObject <- object
  newObject$x <- newX
  newObject$xbar <- matrix(apply(newObject$x[keep,-1],2,mean),
                           nrow=n,ncol=d,byrow=TRUE)
  dimnames(newObject$xbar) <- dimnames(object$xbar)
  
  ## for Beta?
  if(haveBeta){
    newObject$beta <- newBeta
    newObject$betabar <- matrix(apply(newBeta[keep,-1],2,mean),
                                nrow=m,ncol=d+1,byrow=TRUE)
    dimnames(newObject$betabar) <- dimnames(object$betabar)
  }
    
  newObject
}

affineTrans <- function(x,target){
  d <- dim(x)[2]
  x0 <- cbind(x,1)
  if(d>1){
    zeroMat <- 0*x0
    A <- rbind(cbind(x0,zeroMat),
               cbind(zeroMat,x0))
    b <- as.vector(target)
  }
  if(d==1){
    A <- x0
    b <- target
  }
  foo <- solve(A)%*%b
  foo <- matrix(foo,nrow=d+1)
  foo <- cbind(foo,
               c(rep(0,d),1))
  foo
}
