"hurdle" <-
function(count= y ~ .,
                     x = ~1,
                     z = ~1,
                     data=list(),
                     link="logit",
                     dist="poisson",
                     method="BFGS",
                     trace=FALSE,
                     maxit=50000,
                     na.action=na.omit)
{
  cat("Hurdle Count Model\n")
  if(link=="probit"){
      linkfn <- pnorm
      cat("Using probit to model zero vs non-zero\n")
  }
  else{
      linkfn <- plogis
      cat("Using logit to model zero vs non-zero\n")
  }
  
  hurdlePoisson <- function(parms){
      gamma <- parms[1:kz]
      beta <- parms[(kz+1):(kz+kx)]
      
      ## logitfor y=0
      mu <- Z%*%gamma
      p0 <- linkfn(mu)
      
      ## poisson for y>0
      lambda <- exp(X%*%beta)
      f2 <- dpois(Y,lambda=lambda)       ## Poisson prob, actual outcomes
      f2.0 <- dpois(0,lambda=lambda)     ## Poisson prob, zero outcome  
      p1 <- (1-p0)/(1-f2.0) * f2
      
      p0 <- p0[zeroCount]
      p1 <- p1[posCount]
      
      ## log-likelihood
      llh <- sum(log(p0)) + sum(log(p1))
      llh
  }
  
  NegBin <- function(z,lambda,theta){
    arg1 <- lgamma(theta + z) - lgamma(z+1) - lgamma(theta)
    logR <- log(lambda) - log(lambda+theta)
    oneMinusR <- 1 - exp(logR)
    arg2 <- log(oneMinusR)*theta
    arg3 <- logR*z
    out <- arg1 + arg2 + arg3
    out <- exp(out)
    out
  }
  
  hurdleNegBin <- function(parms){
      gamma <- parms[1:kz]
      beta <- parms[(kz+1):(kz+kx)]
      theta <- exp(parms[kz+kx+1])
      
      ## logitfor y=0
      mu <- Z%*%gamma
      p0 <- linkfn(mu)
      
      ## negbin for y>0
      lambda <- exp(X%*%beta)
      
      f2 <- NegBin(Y,lambda,theta)   ## NegBin prob, actual outcomes
      f2.0 <- NegBin(0,lambda,theta) ## NegBin prob, zero outcome  
      p1 <- (1-p0)/(1-f2.0) * f2
      
      p0 <- p0[zeroCount]
      p1 <- p1[posCount]
      
      ## log-likelihood
      llh <- sum(log(p0)) + sum(log(p1))
      llh
  }
  
  hurdleNegBinNoExp <- function(parms){
      gamma <- parms[1:kz]
      beta <- parms[(kz+1):(kz+kx)]
      theta <- parms[kz+kx+1]
      
      ## logitfor y=0
      mu <- Z%*%gamma
      p0 <- linkfn(mu)
      
      ## negbin for y>0
      lambda <- exp(X%*%beta)
      
      f2 <- NegBin(Y,lambda,theta)   ## NegBin prob, actual outcomes
      f2.0 <- NegBin(0,lambda,theta) ## NegBin prob, zero outcome  
      p1 <- (1-p0)/(1-f2.0) * f2
      
      p0 <- p0[zeroCount]
      p1 <- p1[posCount]
      
      ## log-likelihood
      llh <- sum(log(p0)) + sum(log(p1))
      llh
  }
  
  if(dist=="poisson"){
      cat("Using Poisson for counts\n")
      llhfunc <- hurdlePoisson
  }
  else if(dist=="negbin"){
      cat("Using Negative Binomial for counts\n")
      llhfunc <- hurdleNegBin
  }
  else
      stop(paste("distribution",dist,"unsupported"))

  if (method=="Nelder-Mead")
      control <- list(maxit=maxit,trace=trace,REPORT=1,fnscale=-1)
  else if(method=="BFGS")
      control <- list(trace=trace,REPORT=1,maxit=maxit,fnscale=-1,rel.tol=1e-24)
  else
      stop(paste("optimization method",method,"unsupported"))
  
  ## set-up, largely borrowed from glm
  cl <- match.call()
    
  if (missing(data)) 
      data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  mX <- match(c("x", "data", "subset", "weights", "na.action", "offset"),
              names(mf), 0)
  mZ <- match(c("z", "data", "subset", "weights", "na.action", "offset"),
              names(mf), 0)
  mY <- match(c("count", "data", "subset", "weights", "na.action", "offset"),
              names(mf), 0)
  
  mfX <- mf[c(1,mX)]
  names(mfX)[names(mfX)=="x"] <- "formula"
  mfX$drop.unused.levels <- TRUE
  mfX[[1]] <- as.name("model.frame")
  mfX <- eval(mfX, parent.frame())
  mtX <- attr(mfX, "terms")
  
  mfZ <- mf[c(1,mZ)]
  names(mfZ)[names(mfZ)=="z"] <- "formula"
  mfZ$drop.unused.levels <- TRUE
  mfZ[[1]] <- as.name("model.frame")
  mfZ <- eval(mfZ, parent.frame())
  mtZ <- attr(mfZ, "terms")
  
  mfY <- mf[c(1, mY)]
  names(mfY)[names(mfY)=="count"] <- "formula"
  mfY$drop.unused.levels <- TRUE
  mfY[[1]] <- as.name("model.frame")
  mfY <- eval(mfY, parent.frame())
  Y <- model.response(mfY, "numeric")
  if (length(dim(Y)) == 1) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
          names(Y) <- nm
  }
  if(is.null(Y) | length(Y)==0)
      stop("No observations in y")
  
  n <- length(Y)
  tab <- table(Y,exclude=NULL)
  cat("dependent variable y:\n")
  print(tab)
  y0 <- sort(unique(Y))           ## unique values of y
  if(any(floor(y0)!=y0))
      stop("Non-integer values of y encountered, invalid for zeroinfl/count model\n")
  if(min(y0)!=0)
      stop("Minimum value of y is not zero, invalid for zeroinfl/count model\n")
  
  zeroCount <- Y==0
  posCount <- Y>0
  
  X <- if(!is.empty.model(mtX))
      model.matrix(mtX,mfX)
  else
      matrix(, NROW(Y), 0)
  kx <- dim(X)[2]
  
  Z <- if(!is.empty.model(mtZ))
      model.matrix(mtZ,mfZ)
  else
      matrix(, NROW(Y), 0)
  kz <- dim(Z)[2]
  

  ## start values
  cat("generating start values...")
  model0 <- glm.fit(Z,1-pmin(Y,1),family=binomial())
  g0 <- model0$coef
  model1 <- glm.fit(X,Y,family=poisson())
  g1 <- model1$coef
  stval <- c(g0,g1)
  if(dist=="negbin")
      stval <- c(stval,0)
  cat("done\n")

  options(warn=0)                ## suppress warnings
  if(dist=="poisson"){
      foo <- optim(fn=llhfunc,
                   par=stval,
                   method=method,
                   control=control,
                   hessian=TRUE)
  }
  if(dist=="negbin"){
      foo <- optim(fn=llhfunc,
                   par=stval,
                   method=method,
                   control=control,
                   hessian=FALSE)
      stval <- foo$par
      stval[length(stval)] <- exp(stval[length(stval)])
      foo <- optim(fn=hurdleNegBinNoExp,
                   par=stval,
                   method=method,
                   control=control,
                   hessian=TRUE)
  }
  if(foo$convergence != 0)
      stop("optimization failed to converge")
  
  out <- list()
  out$stval <- stval
  out$par <- foo$par
  if(dist=="poisson")
    names(out$par) <- c(dimnames(Z)[[2]],dimnames(X)[[2]])
  if(dist=="negbin")
    names(out$par) <- c(dimnames(Z)[[2]],dimnames(X)[[2]],"theta")
  out$hessian <- foo$hessian
  out$call <- cl
  out$method <- method
  out$dist <- dist
  out$llh <- foo$value
  out$y <- Y
  out$x <- X
  out$z <- Z

  ## need levels and terms to make predict work properly
  out$xlevels <- .getXlevels(mtX,mfX)
  out$zlevels <- .getXlevels(mtZ,mfZ)
  out$xTerms <- mtX
  out$zTerms <- mtZ
  out$link <- link
  out$linkfn <- linkfn
  out$n <- n
  out$kx <- kx
  out$kz <- kz

  class(out) <- c("hurdle")
  out
}

