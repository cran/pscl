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
                  normalize=FALSE,
                  meanzero=normalize,
                  priors=NULL,
                  startvals="eigen",
                  store.item=FALSE,
                  file=NULL,
                  verbose=FALSE){

  cat("ideal: analysis of roll call data via Markov chain Monte Carlo methods.\n\n")

  ## name of rollcall object as unevaluated, parsed expression
  cl <- match.call()
  if(is.null(cl$d))
    cl$d <- d
  if(is.null(cl$codes))
    cl$codes <- codes
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
  if(is.null(cl$normalize))
    cl$normalize <- normalize

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

  if(!is.null(normalize) & d>1){
    cat("normalize option is only meaningful when d=1\n")
  }

  if(normalize != meanzero){
    normalize <- meanzero
    cat("meanzero option is being phased out; normalize provides the same functionality\n")
    cat(paste("For now, we will use your supplied value of meanzero, proceeding with normalize=",
              meanzero,"\n"))
  }
              
  ## pre-process rollcall object
  tmpObject <- object
  if(!is.null(codes)){
    tmpObject$codes <- codes
    if(checkCodes(tmpObject$codes))
      stop("supplied codes fail redundancy checks")
  }
  if(!is.null(dropList)){
    if(verbose)
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
  if(verbose){
    printCodes(codes)
    cat("\n")
  }
  
  if(checkVotes(y$votes,codes))
    stop("rollcall: can't map all votes using supplied codes")

  v <- convertCodes(y,codes)  ## convert to zeros and ones and NAs
  
  ## using a file for storage
  usefile <- !is.null(file)      

  ## check to see how much information will need to be stored
  numrec <- (maxiter-burnin)/thin+1
  if (interactive() &
      ((store.item)&&((n+m)*d*numrec>2000000))
      ||
      ((!store.item)&&((n*d*numrec)>2000000))
      ){
    ans <- readline(paste("The current call to ideal will result in a large object that\n",
                          "will take up a large amount of memory.  Do you want to\n",
                          "continue with the current configuation? (y/n): ",
                          sep=""))
    
    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated execution of ideal.")
  }

  if (interactive() & numrec>1000) {
    ans <- readline(paste("You are attempting to save ",numrec," iterations.  This\n",
                          "could result in a very large object and cause memory problems.\n",
                         "Do you want to continue with the current call to ideal? (y/n): ",
                         sep=""))
    
    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated execution of ideal.")
  }

  cat(paste("Ideal Point Estimation\n\nNumber of Legislators\t\t",
              n,"\nNumber of Items\t\t\t", m, "\n\n"))
  xp <- xpv <- bp <- bpv <- NULL

  ####################################################################
  ## check priors
  ####################################################################
  if(verbose)
    cat("checking for any user-supplied priors...\n")
  if(!is.null(priors)){
    if(!is.list(priors))
      stop("priors must be a list")
    
    if(all(unlist(lapply(priors,is.null))))
      stop("priors supplied in a list, but all elements are NULL")
    
    if(sum(unlist(lapply(priors,is.na)))>0)
      stop("priors contain missing values, which is not allowed")
  
    ## now check individual elements of prior list
    if(!is.null(priors$xp)){
      if(length(priors$xp)==1)     ## user supplied a scalar
        xp <- matrix(priors$xp,n,d)
      ## coerce a vector to a matrix
      if(length(priors$xp)>1 & d==1 & !is.matrix(priors$xp))
        xp <- matrix(priors$xp,n,d)
      if(is.matrix(priors$xp))
        xp <- priors$xp
    }
    else{
      if(verbose)
        cat("no prior means supplied for ideal points,\n",
            "setting to default of 0\n")
      xp <- matrix(0,n,d)
    }
    
    if(!is.null(priors$xpv)){
      if(length(priors$xpv)==1)    ## user supplied a scalar
        xpv <- matrix(priors$xpv,n,d)
      ## coerce a vector to a matrix
      if(length(priors$xpv)>1 & d==1 & !is.matrix(priors$xpv))
        xpv <- matrix(priors$xpv,n,d)
      if(is.matrix(priors$xpv))
        xpv <- priors$xpv
    }
    else{
      if(verbose)
        cat("no prior precisions supplied for ideal points,\n",
            "setting to default of 1\n")
      xpv <- matrix(1,n,d)
    }
    
    if(!is.null(priors$bp)){
      if(length(priors$bp)==1)    ## user supplied a scalar
        bp <- matrix(priors$bp,m,d+1)
      if(is.matrix(priors$bp))
        bp <- priors$bp
    }
    else{
      if(verbose)
        cat("no prior means supplied for item parameters,\n",
            "setting to default to 0\n")
      bp <- matrix(0,m,d+1)
    }
    
    if(!is.null(priors$bpv)){
      if(length(priors$bpv)==1)   ## user supplied a scalar
        bpv <- matrix(priors$bpv,m,d+1)
      if(is.matrix(priors$bpv))
        bpv <- priors$bpv
    }
    else{
      if(verbose)
        cat("no prior precisions supplied for item parameters,\n",
            "setting to default of .01\n")
      bpv <- matrix(.01,m,d+1)
    }
    
    if (((nrow(xp) != n)||(ncol(xp) != d)) || ((nrow(xpv)!=n)||(ncol(xpv)!=d))) {
      stop("Dimensions of xp or xpv not n by d")
    }
    
    if (((nrow(bp) != m)||(ncol(bp) != (d+1))) || ((nrow(bpv)!=m)||(ncol(bpv)!=(d+1)))) {
      stop("Dimensions of bp or bpv not m by d+1")
    }
  }
  
  ## ##################################################################
  ## if we get this far with priors still NULL
  ## then revert to defaults
  ## ##################################################################
  if(is.null(xp)){
    if(verbose)
      cat("setting prior means for ideal points to all zeros\n")
    xp <- matrix(0,n,d)
  }
  if(is.null(xpv)){
    if(verbose)
      cat("setting prior precisions for ideal points to all 1\n")
    xpv <- matrix(1,n,d)
  }
  if(is.null(bp)){
    if(verbose)
      cat("setting prior means for item parameters to all zeros\n")
    bp <- matrix(0,m,d+1)
  }
  if(is.null(bpv)){
    if(verbose)
      cat("setting prior precisions for item parameters to all 0.01\n")
    bpv <- matrix(0.01,m,d+1)
  }
  
  xp <- as.vector(t(xp))
  xpv <- as.vector(t(xpv))
  bp <- as.vector(t(bp))
  bpv <- as.vector(t(bpv)) 

  ################################################################
  ## check for start values - create if not supplied
  ################################################################
  if(verbose)
    cat("\nchecking start values...\n")
  xstart <- bstart <- NULL
  options(warn=-1)

  if(!is.list(startvals)){
    if(startvals=="eigen" | is.null(startvals)){
      xstart <- x.startvalues(v,d=d,verbose=verbose)
      bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
      bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                       bstart,
                       bp + 2*sign(bstart-bp)/sqrt(bpv))
    }

    if(startvals=="random"){
      if(verbose)
        cat("generating start values for ideal points by iid sampling from N(0,1)\n")
      xstart <- matrix(rnorm(n*d),n,d)
      bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
      bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                       bstart,
                       bp + 2*sign(bstart-bp)/sqrt(bpv))
    }
  }

  ## user has passed something in startvals
  if(is.list(startvals)){
    if(!is.null(startvals$xstart)){
      if(length(startvals$xstart) != n*d)
        stop("length of xstart not n by d")
      if(d==1)
        xstart <- matrix(startvals$xstart,ncol=1)
      else
        xstart <- startvals$xstart
      if (sum(is.na(xstart))!=0)
        stop("xstart contains missing values")
    }
    
    if(!is.null(startvals$bstart)){
      if(length(startvals$bstart) != m*(d+1))
        stop("length of bstart not m by d+1")
      bstart <- startvals$bstart
      if(sum(is.na(bstart))!=0)
        stop("bstart contains missing values")
    } 
  } 

  ## final check
  if(is.null(xstart)){
    cat("no user-supplied start values found\n")
    xstart <- x.startvalues(v,d,verbose=TRUE)
  }
  if(is.null(bstart)){
    bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
    bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                     bstart,
                     bp + 2*sign(bstart-bp)/sqrt(bpv))
  }

  ## report to user
  if(verbose){
    cat("using the following start values for ideal points (summary follows):\n")
    print(summary(xstart))
    
    cat("using the following start values for item parameters (summary follows):\n")
    print(summary(bstart))
  }

  xstart <- as.vector(t(xstart))
  bstart <- as.vector(t(bstart))
  
  options(warn=0)

  ##############################################################
  ## end error checking
  ##############################################################

  yToC <- ifelse(is.na(v), 9, v)
  yToC <- as.vector(t(yToC))
  cat("\nStarting MCMC Iterations...\n")

  ## ############################################
  ## two versions, one with usefile option
  ## ############################################
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
      cat(",",
          paste("\"",
                c(paste("b",
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
                 as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=NULL,
                 boutput=NULL,as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file),
                 as.integer(verbose))
  }
  ## not saving output to file, saving output to memory
  else if (!store.item) {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC), 
                 as.integer(maxiter), as.integer(thin), as.integer(impute),      
                 as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=NULL,as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file),
                 as.integer(verbose))
  }
  else {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC),          
                 as.integer(maxiter), as.integer(thin), as.integer(impute),
                 as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=as.double(rep(0,m*(d+1)*numrec)),as.integer(burnin),
                 as.integer(usefile), as.integer(store.item), as.character(file),
                 as.integer(verbose))
  }

  cat("\n")

  ## parse returns from C job
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
    if(verbose)
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
    if(verbose)
    cat("done\n")
  }
  
  if(store.item & !is.null(b)){
    if(verbose)
      cat("and for bill parameters...")
    keep <- b[,1] > burnin
    betabar <- apply(b[keep,-1],2,mean)
    betabar <- matrix(betabar,m,d+1,byrow=TRUE)
    if(verbose)
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

  ## and, finally, if the user wanted meanzero
  if(normalize)
    out <- postProcess(out,
                       constraints="normalize")
  
  return(out) 
}


gencolnames <- function(name, d, beta=FALSE) {
  if(d>1){            ## more than one dimension?
    dname <- NULL
    for(i in 1:d){
      dname <- c(dname,paste(name,"d",i,sep=""))
    }
    if(beta)
      dname <- c(dname,paste(name,"Difficulty",sep=""))
    dname <- matrix(dname,ncol=length(name),byrow=TRUE)
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

x.startvalues <- function(x,d,scale=TRUE,constraint=NULL,verbose=FALSE){
  if(verbose)
    cat("will use eigen-decomposition method to get start values for ideal points...")

  ## from Jong Hee Park
  row.mean <- apply(x, 1, mean, na.rm=TRUE)
  col.mean <- apply(x, 2, mean, na.rm=TRUE)
  dc1 <- sweep(x, 1, row.mean)
  dc2 <- sweep(dc1, 2, col.mean)
  dc <- dc2 + mean(x, na.rm = T)
  
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
  if(verbose)
    cat("done\n")
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

b.startvalues <- function(v,x,d,verbose=FALSE){
  m <- dim(v)[2]
  if(verbose)
    cat(paste("running",
              m,
              "vote-specific probit GLMs\n",
              "for start values for item/bill parameters\n",
              "conditional on start values for ideal points..."))
  
  b <- matrix(NA,m,d+1)
  for(j in 1:m){
    b[j,] <- probit(y=v[,j],x=x)
  }
  b[,d+1] <- -b[,d+1]     ## flip the sign on the intercept, make it a difficulty parameter
  if(verbose)
    cat("done\n")
  b
}

