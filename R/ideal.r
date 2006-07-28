## IDEAL

## ## for creating a package
.package.Name <- "pscl"

## .First.lib<-function(lib,pkg) {
##   library.dynam(.package.Name,
##                 pkg,
##                 lib)
## }

ideal <- function(object,
                  codes=object$codes,
                  dropList=list(codes="notInLegis",lop=0),
                  d=1,
                  maxiter=10000,
                  thin=100,
                  burnin=5000,
                  impute=FALSE,
                  meanzero=FALSE,
                  priors=NULL,
                  startvals=NULL,
                  store.item=FALSE,
                  file=NULL){

  cat("ideal: analysis of roll call data via Markov chain Monte Carlo methods.\n\n")

  ## name of rollcall object as unevaluated, parsed expression
  cl <- match.call()
  if(is.null(cl$d))
    cl$d <- d
  if(is.null(cl$codes))
    cl$d <- codes
  if(is.null(cl$dropList))
    cl$dropList <- dropList
  if(is.null(cl$maxiter))
    cl$maxiter <- maxiter
  if(is.null(cl$thin))
    cl$thin <- thin
  if(is.null(cl$burnin))
    cl$burnin <- burnin
  if(is.null(cl$impute))
    cl$impute <- impute
  if(is.null(cl$store.item))
    cl$store.item <- store.item
  if(is.null(cl$meanzero))
    cl$meanzero <- meanzero

  localMeanZero <- meanzero   ## copy user-supplied

  ## check validity of user arguments
  if (!("rollcall" %in% class(object)))
    stop("object must be of class rollcall")
  if(((d%%1) != 0) || (d<1)){
    stop("d is not a positive integer")
  }

  if(((thin%%1)!=0) || (thin<1)) {
    stop("thin is not a positive integer")
  }

  if(((maxiter%%1)!=0) || (maxiter<1)) {
    stop("maxiter is not a positive integer")
  }

  if(!is.list(dropList))
    stop("dropList must be a list")

  if(!is.list(codes))
    stop("codes must be a list")

  ##check iterations and thinning
  if ((maxiter%%thin)!=0) {
    stop("maxiter must be a multiple of thin")
  }

  if ((burnin%%thin)!=0) {
    stop("burnin must be a multiple of thin")
  }

  if (burnin >= maxiter)
    stop("burnin must be less than maxiter")

  
  ## pre-process rollcall object
  tmpObject <- object
  if(!is.null(codes)){
    tmpObject$codes <- codes
    if(checkCodes(tmpObject$codes))
      stop("supplied codes fail redundancy checks")
  }
  if(!is.null(dropList)){
    cat(paste("Subsetting rollcall object",
              as.name(cl$object),
              "using dropList\n"))
    y <- dropRollCall(tmpObject,dropList)  ## any subsetting to do?
  }
  else
    y <- tmpObject
  rm(tmpObject)
  
  n <- dim(y$votes)[1]
  m <- dim(y$votes)[2]
  legis.names <- dimnames(y$votes)[[1]]
  vote.names <- dimnames(y$votes)[[2]]

  ## map roll call votes into binary format required by ideal
  printCodes(codes)
  cat("\n")

  if(checkVotes(y$votes,codes))
    stop("rollcall: can't map all votes using supplied codes")

  v <- convertCodes(y,codes)  ## convert to zeros and ones and NAs
  
  ## using a file for storage
  usefile <- !is.null(file)      

  ## check to see how much information will need to be stored
  numrec <- (maxiter-burnin)/thin+1
  if (((store.item)&&((n+m)*d*numrec>2000000))
       ||((!store.item)&&((n*d*numrec)>2000000))){
    ans <- readline(paste("The current call to ideal will result in a large object that ",
                          "will take up a large amount of memory.  Do you want to ",
                          "continue with the current parameter values? (y/n): ",
                          sep=""))

    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated ideal run.")
  }

  if (numrec>1000) {
    ans <- readline(paste("You are attempting to save ",numrec," iterations.  This ",
                         "could result in a very large object and cause memory problems.  ",
                         "Do you want to continue with the current call to ideal? (y/n): ",
                         sep=""))

    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated ideal run.")
  }

  cat(paste("Ideal Point Estimation\n\nNumber of Legislators\t\t",
              n,"\nNumber of Items\t\t\t", m, "\n\n"))
  xp <- xpv <- bp <- bpv <- NULL

  ## check priors
  if(!is.null(priors)){
    if(!is.list(priors))
      stop("priors must be a list with elements xp, xpv, bp, bpv")

    if(sum(is.na(match(c("xp","xpv","bp","bpv"),names(priors))))!=0)
      stop("priors must be a list with elements xp, xpv, bp, bpv")

    localMeanZero <- FALSE
    xp <- as.matrix(priors$xp)
    xpv <- as.matrix(priors$xpv)
    bp <- as.matrix(priors$bp)
    bpv <- as.matrix(priors$bpv)
    if (sum(is.na(c(xp,xpv,bp,bpv)))!=0) {
      stop("Priors contain missing values, which is not allowed")
    }

    if (((nrow(xp) != n)||(ncol(xp) != d)) || ((nrow(xpv)!=n)||(ncol(xpv)!=d))) {
      stop("Dimensions of xp or xpv not n by d")
    }

    if (((nrow(bp) != m)||(ncol(bp) != (d+1))) || ((nrow(bpv)!=m)||(ncol(bpv)!=(d+1)))) {
      stop("Dimensions of bp or bpv not m by d+1")
    }
  }
  else {   ## no user-supplied priors
    xp <- matrix(rep(0, n*d), nrow=n)
    xpv <- matrix(rep(1, n*d), nrow=n)
    bp <- matrix(rep(0,m*(d+1)), nrow=m)
    bpv <- matrix(rep(0.01, m*(d+1)), nrow=m)
  }

  xp <- as.vector(t(xp))
  xpv <- as.vector(t(xpv))
  bp <- as.vector(t(bp))
  bpv <- as.vector(t(bpv)) 

  ################################################################
  ## check for start values - create if not supplied
  ################################################################
  xstart <- bstart <- NULL
  options(warn=-1)
  if (!is.null(startvals)) {
    if(!is.list(startvals))
      stop("startval must be a list with elements xstart and bstart")
    if(sum(is.na(match(c("xstart","bstart"),names(startvals))))!=0)
      stop("startval must be a list with elements xstart and bstart")
    xstart <- as.matrix(startvals$xstart)

    if ((nrow(xstart) != n)||(ncol(xstart) != d)){
      stop("Dimensions of xstart not n by d")
    }

    if (sum(is.na(xstart))!=0)
      stop("xstart contains missing values")

    bstart <- as.matrix(startvals$bstart)

    if ((nrow(bstart) != m)||(ncol(bstart) != (d+1))) {
      stop("Dimensions of bstart not m by d+1")
    }

    if(sum(is.na(bstart))!=0)
      stop("bstart contains missing values")
  } 
  else {
    if((n>500) | (!is.null(priors))){
      if (n>500)
        cat(paste("n =",n,
                  "is too many subjects/legislators for eigen-decomposition\n",
                  "for start values (this option for x.start capped at n=500)\n"))
      cat("setting start values to zero for all legislators\n")
      xstart <- matrix(0,n,d)
      bstart <- matrix(0,m,d+1)
    }
    else
      {
        xstart <- x.startvalues(v,d=d)
        cat(paste("running",m,"vote-specific probits\nfor start values for item/bill parameters\n"))
        bstart <- b.startvalues(v,xstart,d=d)
        bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv), bstart, bp + 2*sign(bstart-bp)/sqrt(bpv))
      }
  }

  xstart <- as.vector(t(xstart))
  bstart <- as.vector(t(bstart))
  options(warn=0)

  ##############################################################
  ## end error checking
  ##############################################################

  yToC <- ifelse(is.na(v), 9, v)
  yToC <- as.vector(t(yToC))
  cat("\nStarting iterations...\n")
  if (usefile) {
    if (length(legis.names) == n) {
      cat(paste("\"",c("Iteration",legis.names),"\"", sep="", collapse=","),
          file=file)
    }
    else {
      cat(paste("\"",c("Iteration",paste("x", 1:n, sep="")),"\"",
                sep="", collapse=","),
          file=file)
    }

    if (store.item){
            cat(",", paste("\"",c(paste("b",
                                        as.vector(apply(expand.grid(1:m,1:(d+1)),1,paste,collapse=".")),
                                        sep=".")),"\"",
                           sep="", collapse=","),
                sep="", file=file, append=TRUE)
          }
    cat("\n", file=file, append=TRUE)
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC), 
                 as.integer(maxiter), as.integer(thin), as.integer(impute),
                 as.integer(localMeanZero), as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=NULL,
                 boutput=NULL,as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file))
  }

  else if (!store.item) {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC), 
                 as.integer(maxiter), as.integer(thin), as.integer(impute),      
                 as.integer(localMeanZero), as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=NULL,as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file))
  }
  else {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC),          
                 as.integer(maxiter), as.integer(thin), as.integer(impute),
                 as.integer(localMeanZero), as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=as.double(rep(0,m*(d+1)*numrec)),as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file))
  }

  cat("\n")

  ## read output
  if (!usefile) {
    x <- output$xoutput
    x <- matrix(x,nrow=numrec,byrow=T)
    itervec <- seq(burnin,maxiter,by=thin)
    x <- cbind(itervec,x)
    rownames(x) <- x[,1]
    colnames(x) <- gencolnames(legis.names,d)
    if (store.item) {
      b <- output$boutput
      b <- matrix(b,nrow=numrec,byrow=T)
      b <- cbind(itervec,b)
      rownames(b) <- b[,1]
      colnames(b) <- gencolnames(vote.names,d,beta=T)
    }
    else {
      b <- NULL
    }
  }
  else {   ## output went to a file
    b <- x <- NULL
  }

  ## compute some summary stats now, since we almost always use them
  xbar <- betabar <- NULL
  if(!is.null(x)){
    cat("MCMC sampling done, computing posterior means for ideal points...\n")
    keep <- x[,1] > burnin
    xbar <- apply(x[keep,-1],2,mean)
    xbar <- matrix(xbar,n,d,byrow=TRUE)
    mnames <- NULL
    if(d>1){
      for(j in 1:d)
        mnames <- c(mnames,paste("Dimension",j))
    }
    dimnames(xbar) <- list(legis.names,mnames)
    cat("done\n")
  }
  
  if(store.item){
    cat("and for bill parameters...")
    betabar <- apply(b[keep,-1],2,mean)
    betabar <- matrix(betabar,m,d+1,byrow=T)
    cat("done\n")
  }

  ## wrap up for return to user
  out <- list(n=n,m=m,d=d,
              codes=codes,
              x=x,
              beta=b,
              xbar=xbar,
              betabar=betabar,
              call=cl)

  class(out) <- c("ideal")
  return(out) 
}


gencolnames <- function(name, d, beta=F) {
  if(d>1){            ## more than one dimension?
    dname <- NULL
    for(i in 1:d){
      dname <- c(dname,paste(name,"d",i,sep=""))
    }
    if(beta)
      dname <- c(dname,paste(name,"Intercept",sep=""))
    dname <- matrix(dname,ncol=length(name),byrow=T)
    dname <- as.vector(dname)
    dname <- c("Iteration",dname)
  }
  else {
    if(beta){
      if(beta)
        dname <- c(name,paste(name,"Intercept",sep=""))
      dname <- matrix(dname,ncol=length(name),byrow=T)
      dname <- as.vector(dname)
      dname <- c("Iteration",dname)
    }
    else {
      dname <- c("Iteration",name)
    }
  }
  dname
}

x.startvalues <- function(x,d,scale=TRUE,constraint=NULL){
  dc <- apply(x,2,function(x)x-mean(x,na.rm=T))
  dc <- apply(dc,2,function(x)x-mean(x,na.rm=T))
  dc <- dc + mean(x,na.rm=T)
  r <- cor(t(dc),use="pairwise")
  r[is.na(r)] <- 0
  e <- eigen(r)
  v <- e$vectors[,1:d]
  v <- as.matrix(v)
  if(scale){
    for(i in 1:d){
      v[,i] <- v[,i]*sqrt(e$value[i])
    }
  }

  if (!is.null(constraint)) {
    v <- predict(lm(constraint ~ v), newdata=as.data.frame(v)) 
  }
  v
}

probit <- function(y,x){
  glmobj <- glm(y ~ x,
                family=binomial(link=probit))
  b <- coef(glmobj)
  k <- length(b)
  b <- b[c(2:k,1)]
  b
}

b.startvalues <- function(v,x,d){
  m <- dim(v)[2]
  b <- matrix(NA,m,d+1)
  for(j in 1:m){
    b[j,] <- probit(y=v[,j],x=x)
  }
  b
}

