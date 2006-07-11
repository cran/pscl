"zeroinfl" <-
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
    cat("Zero-Inflated Count Model\n")
    if(link=="probit"){
        linkfn <- pnorm
        cat("Using probit to model zero vs non-zero\n")
    }
    else{
        linkfn <- plogis
        cat("Using logit to model zero vs non-zero\n")
    }
    
    zeroinflPoisson <- function(parms){
        gamma <- parms[1:kz]
        beta <- parms[(kz+1):(kz+kx)]
        
        phi <- linkfn(Z%*%gamma)
        xb <- X%*%beta
        lambda <- exp(xb)
        
        arg1 <- log(phi + (1-phi)*(exp(-lambda)))
        arg1 <- arg1[zeroCount]
        
        arg2 <- log(1 - phi) - lambda + Y*xb - lgamma(Y+1)
        arg2 <- arg2[posCount]
        
        llh <- sum(arg1) + sum(arg2)
        llh
    }
    
    NegBin <- function(Y,lambda,theta){
                                        #dnbinom(y,mu=lambda,size=theta)
        arg1 <- gamma(theta + Y)/(gamma(Y+1)*gamma(theta))
        r <- lambda/(lambda+theta)
        arg2 <- (1-r)^theta
        arg3 <- r^Y
        out <- arg1*arg2*arg3
        out
    }
    
    LogNegBin <- function(Y,lambda,theta){
                                        #dnbinom(x=y,size=theta,mu=lambda,log=T)
        arg1 <- lgamma(Y+theta)-lgamma(Y+1)-lgamma(theta)
        arg2 <- theta*(log(theta)-log(theta+lambda))
        arg3 <- Y*(log(lambda)-log(theta+lambda))
        out <- arg1+arg2+arg3
        out
    }
    
    zeroinflNegBin <- function(parms){
        gamma <- parms[1:kz]
        beta <- parms[(kz+1):(kz+kx)]
        theta <- exp(parms[(kz+kx)+1])
        
        mu <- Z%*%gamma
        lambda <- exp(X%*%beta)
        
        ## binary link for p(y=0)
        f0 <- linkfn(mu)
        arg2 <- log(1-f0) + theta*(log(theta)-log(theta+lambda))
        arg2 <- exp(arg2)
        ##p0 <- f0 + (1-f0)*(theta/(theta+lambda))^theta
        p0 <- f0 + arg2
        p0 <- p0[zeroCount]
        
        ## negbin for y>0
        ##p1 <- (1-f0)*NegBin(y,lambda,theta)
        ##p1 <- p1[y>0]
        lp1 <- log(1-f0) + LogNegBin(Y,lambda,theta)
        lp1 <- lp1[posCount]
        
        ## log-likelihood
                                        #llh <- sum(log(p0)) + sum(log(p1))
        llh <- sum(log(p0)) + sum(lp1)
        llh
    }
    
    if(dist=="poisson"){
        cat("Using Poisson for counts\n")
        llhfunc <- zeroinflPoisson
    }
    if(dist=="negbin"){
        cat("Using Negative Binomial for counts\n")
        llhfunc <- zeroinflNegBin
    }
    
    if (method=="Nelder-Mead")
        control <- list(maxit=maxit,rel.tol=1e-12,trace=trace,REPORT=1,fnscale=-1)
    else if(method=="BFGS")
        control <- list(trace=trace,REPORT=1,maxit=maxit,fnscale=-1,rel.tol=1e-24)
    else
        stop("invalid method")


    ## set-up, largely borrowed from glm
    cl <- match.call()

    if (missing(data)) 
        data <- environment(count)   ## error spotted by Bettina Gruen <gruen@ci.tuwien.ac.at>
    
    mf <- match.call(expand.dots = FALSE)
    if(is.null(mf$x))
      mf$x <- ~ 1  ## error spotted by Bettina Gruen <gruen@ci.tuwien.ac.at>
    if(is.null(mf$z))
      mf$z <- ~ 1  ## error spotted by Bettina Gruen <gruen@ci.tuwien.ac.at>
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
    m <- length(y0)                 ## number of unique values
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
    
    cat("MLE begins...\n")
    foo <- optim(fn=llhfunc,
                 par=stval,
                 method=method,
                 control=control,
                 hessian=TRUE)
    if(foo$convergence != 0)
        stop("optimization failed to converge")
    cat("done\n")

    out <- list()
    out$stval <- stval
    out$par <- foo$par
    out$theta <- NULL
    if(dist=="poisson")
        names(out$par) <- c(dimnames(Z)[[2]],dimnames(X)[[2]])
    if(dist=="negbin"){
        names(out$par) <- c(dimnames(Z)[[2]],dimnames(X)[[2]],"log(theta)")
        out$theta <- exp(out$par["log(theta)"])
    }
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

    class(out) <- "zeroinfl"
    out
}

