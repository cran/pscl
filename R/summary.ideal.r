printHeaderIdeal <- function(x){
  cat(paste("ideal was called as follows:\n"))
  print(x$call)
  cat("\n")

  cat(paste("Number of Legislators:\t",x$n,"\n"))
  cat(paste("Number of Votes:\t",x$m,"\n"))
  cat(paste("Number of Dimensions:\t",x$d,"\n"))
  cat(paste("Number of Iterations:\t",x$call$maxiter,"\n"))
  cat(paste("\tThinned By:\t",x$call$thin,"\n"))
  cat(paste("\tBurn-in:\t",x$call$burnin,"\n\n"))

  invisible(NULL)
}

## summary and print functions
print.ideal <- function(x, ...) {
  if(class(x) != "ideal")
    stop("object passed to print.ideal is not of class ideal\n")

  cat("Markov chain Monte Carlo Analysis of Roll Call Data\n")
  if(x$d==1)
    cat("       (2-parameter item-response modeling)        \n")
  else
    cat("     (multidimensional item-response modeling)     \n")
  cat("===================================================\n\n")

  printHeader(eval(x$call$object))
  printHeaderIdeal(x)
  
  cat("Ideal Points: Posterior Mean\n")
  print(round(x$xbar,2))
  cat("\n")
  invisible(NULL)
}

summary.ideal <- function(object,
                          quantiles=c(.025,.975),
                          burnin=NULL,
                          sort=TRUE,
                          include.beta=FALSE,
                          ...){

  if(class(object)!="ideal")
    stop("summary.ideal only defined for objects of class ideal")

  if(is.null(burnin))
    keep <- checkBurnIn(object,eval(object$call$burnin))
  else
    keep <- checkBurnIn(object,burnin)
  
  xm <- bm <- xQuantTab <- bQuantTab <- NULL
  xResults <- list()
  bResults <- list()
  bSig <- list()

  ## get quantiles of x
  if(!is.null(quantiles)){
    localQuantiles <- sort(quantiles)
    xQuantTab <- as.matrix(apply(object$x[keep,-1],2,
                                 quantile,
                                 probs=localQuantiles,
                                 na.rm=T))
  }
  if (length(quantiles) > 1)
    xQuantTab <- t(xQuantTab)
  rownames(xQuantTab) <- colnames(object$x)[-1]
  if(is.null(burnin))
    xm <- object$xbar
  else{
    xm <- apply(object$x[keep,-1],2,mean)
    xm <- matrix(xm,ncol=object$d,byrow=TRUE)
  }
  xsd <- apply(object$x[keep,-1],2,sd)
  xQuantTab <- cbind(xsd,xQuantTab)

  ## loop over dimensions
  for(j in 1:object$d){
    thisDimension <- seq(from=j,by=object$d,length=object$n)
    xResults[[j]] <- cbind(xm[,j],
                           xQuantTab[thisDimension,])
    cNames <- c("Mean","Std.Dev.",
                paste(as.character(localQuantiles*100),
                      "%",
                      sep=""))
    dimnames(xResults[[j]])[[2]] <- cNames
    if(sort)
      xResults[[j]] <- xResults[[j]][order(xResults[[j]][,1]),]
  }

  ##################################################################
  ## get beta summaries
  if ((!is.null(object$beta)) && (include.beta)){
    tmpRollCall <- computeMargins(eval(object$call$object),
                                  dropList=object$call$dropList)
    if(!is.null(quantiles))
      bQuantTab <- as.matrix(apply(object$beta[keep,-1],2,
                                quantile,
                                probs=localQuantiles,na.rm=T))
    if (length(quantiles)>1)
      bQuantTab <- t(bQuantTab)
    rownames(bQuantTab) <- colnames(object$beta)[-1]
    if(is.null(burnin))
      bm <- object$betabar
    else{
      bm <- apply(object$beta[keep,-1],2,mean)
      bm <- matrix(bm,ncol=object$d+1,byrow=TRUE)
    }
    bsd <- apply(object$beta[keep,-1],2,sd)
    bQuantTab <- cbind(bsd,bQuantTab)

    ## loop over dimensions
    for(b in 1:(object$d+1)){
      thisDimension <- seq(from=b,by=object$d+1,length=object$m)
      bResults[[b]] <- cbind(bm[,b],
                             bQuantTab[thisDimension,])
      dimnames(bResults[[b]])[[2]][1] <- "Mean"
      dimnames(bResults[[b]])[[2]][2] <- "sd"

      ## if available, tack on margins data
      if(!is.null(tmpRollCall$voteMargins)){
        bResults[[b]] <- cbind(bResults[[b]],
                               tmpRollCall$voteMargins)
      }
    }
    names(bResults) <- c(paste("Discrimination Parameter Dimension",
                               1:object$d),
                         "Intercept")
    
    ## "significance tests" for discrimination parameters
    nQ <- length(quantiles)
    if(nQ==2){
      ## we have confidence interval
      width <- abs(localQuantiles[2] - localQuantiles[1])
      sigFunc <- function(x){
        out <- sign(x[,3])==sign(x[,4])
        labs <- rep(paste(round(width*100),"% CI",sep=""),2)
        labs[1] <- paste(labs[1],"does NOT overlap 0")
        labs[2] <- paste(labs[2],"overlaps 0")
        out <- factor(out,
                      levels=c("TRUE","FALSE"),
                      labels=labs)
        out
      }
      bSig <- lapply(bResults[-length(bResults)],sigFunc)
      names(bSig) <- paste("Discrimination Parameter Dim",1:object$d)
    }
  }

  #####################################################################
  ## summarize by party
  pall.final <- NULL
  party <- eval(object$call$object)$legis.data$partyName
  if(is.null(party))
    party <- eval(object$call$object)$legis.data$party
  if(!is.null(party)){                       ## we have some party info
    nms <- NULL
    for (b in 1:object$d){       ## loop over dimensions
      pall <- NULL
      pm <- tapply(xm[,b],party,mean)
      if (!is.null(quantiles)){  ## quantiles, if requested
        pq <- tapply(xm[,b],party,quantile,probs=localQuantiles)
        for(j in 1:length(pq)){
          pall <- rbind(pall,pq[[j]])
        }
      }
      pall <- cbind(pm,pall)
      pall.final <- rbind(pall.final,pall)
      nms <- c(nms,paste(rownames(pall),": Dimension ",b,sep=""))
    }

    colnames(pall.final)[1] <- "Mean"
    rownames(pall.final) <- nms
  }

  #####################################################################
  ## gather for output
  out <- list(object=match.call()$object,
              xResults=xResults,
              bResults=bResults,
              bSig=bSig,
              party.quant=pall.final,
              sort=sort)

  class(out) <- "summary.ideal"

  out
}

print.summary.ideal <- function(x, digits=3, ...){ 
  if (!("summary.ideal" %in% class(x)))
    stop("object passed to print.summary.ideal must be of class summary.ideal")

  cat("Markov chain Monte Carlo Analysis of Roll Call Data\n")
  m <- eval(x$object)$m
  d <- eval(x$object)$d
  if(d==1)
    cat("       (2-parameter item-response modeling)        \n")
  else
    cat("     (multidimensional item-response modeling)     \n")

  printHeader(eval(eval(x$object)$call$object))
  printHeaderIdeal(eval(x$object))
  
  if(!is.null(x$party.quant)) {
    cat("Ideal Points (Posterior Means), by Party\n")
    print(round(x$party.quant,digits))
    cat("\n")
  }  

  for(j in 1:d){
    if(x$sort)
      cat(paste("Ideal Points, Dimension ",j,
                "(sorted by posterior means):\n",sep=""))
    else
      cat(paste("Ideal Points, Dimension ",j,":\n",sep=""))
    print(round(x$xResults[[j]],digits))
    cat("\n")
  }

  ## report statistical tests of significance
  if(length(x$bSig)!=0){
    cat("Statistical tests of discrimination parameters:\n")
    if(d==2){ ## do a cross-tabulation
      cat("dimension 1 (rows) against dimension 2 (columns)\n")
      print(table(x$bSig[[1]],x$bSig[[2]]))
    }
    else{
      for(j in 1:d){
        cat("Dimension:",j)
        print(table(x$bSig[[j]]))
        cat("\n")
      }
    }
  }

  if(length(x$bResults)!=0){
    for(j in 1:d){
      cat(paste(names(x$bResults)[j],":\n"))
      theseResults <- x$bResults[[j]]
      foo <- x$bSig[[j]] == (levels(x$bSig[[j]])[2])
      fooChar <- rep("  ",m)
      fooChar[foo] <- "NS"
      dimnames(theseResults)[[1]] <- paste(dimnames(theseResults)[[1]],
                                           fooChar)
      print(round(theseResults,digits))
      cat("\n")
    }
    cat(paste(names(x$bResults)[d+1],":\n"))
    print(round(x$bResults[[d+1]],digits))
  }

  cat("\n")
  invisible(NULL)
}







