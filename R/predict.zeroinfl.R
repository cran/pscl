"predict.zeroinfl" <-
function(object,
         newdata,
         se.fit=FALSE,conf=.95,MC=1000,
         type=c("response","prob"),
         na.action=na.pass,
         ...)
{
    if(class(object)!="zeroinfl")
        stop("error: predict.zeroinfl only takes objects of class zeroinfl\n")


    localType <- match.arg(type)
    linkfn <- object$linkfn
    parms <- object$par
    kx <- object$kx
    kz <- object$kz
    dist <- object$dist
    gamma <- parms[1:kz]
    beta <- parms[(kz+1):(kz+kx)]
    theta <- object$theta
    out <- list()

    ## if no new data supplied
    if(missing(newdata)){
        X <- object$x
        Z <- object$z
    }
    else{
        mX <- model.frame(object$xTerms,newdata,na.action=na.action,
                          xlev = object$xlevels)
        mX <- model.matrix(object$xTerms,mX,contrasts=attr(object$x,"contrasts"))
        mZ <- model.frame(object$zTerms,newdata,na.action=na.action,
                          xlev = object$zlevels)
        mZ <- model.matrix(object$zTerms,mZ,contrasts=attr(object$z,"contrasts"))
        if(dim(mX)[1]!=dim(mZ)[1])
            stop(paste("Different numbers of rows in newdata generated ",
                       "for count and zero-inflated parts of the model",
                       sep="\n"))
        else{
            X <- mX
            Z <- mZ
        }
     }   

    mu <- exp(X%*%beta)
    phi <- linkfn(Z%*%gamma)
    yhat <- mu - mu*phi
    
    ## if user requests standard errors, use simulation
    if(se.fit){
        loadLib <- try(library(mvtnorm),silent=TRUE)
        if(inherits(loadLib,"try-error"))
            stop("could not load library mvtnorm\n",
                 "standard errors and confidence intervals can't be computed")

        vc <- try(solve(-object$hessian),silent=TRUE)
        if(inherits(vc,"try-error"))
            stop(paste("standard errors and confidence intervals can not be computed",
                       "because hessian matrix of MLEs can not be inverted"),
                 sep="\n")

        yhat.sim <- matrix(NA,MC,dim(X)[1])
        cat("commencing Monte Carlo simulations for predicted counts\n")
        for(i in 1:MC){
            cat(paste("MC iterate",i,"of",MC,"\n"))
            parms.sim <- rmvnorm(n=1,mean=parms,sigma=vc)
            gamma <- parms.sim[1:kz]
            beta <- parms.sim[(kz+1):(kz+kx)]
            mu.sim <- exp(X%*%beta)
            phi.sim <- linkfn(Z%*%gamma)
            yhat.sim[i,] <- mu.sim - mu.sim*phi.sim
        }
        out$lower <- apply(yhat.sim,2,quantile,(1-conf)/2)
        out$upper <- apply(yhat.sim,2,quantile,1-((1-conf)/2))
        out$se <- apply(yhat.sim,2,sd)
    }

    ## user requests predicted probabilities
    if(localType=="prob"){
        yUnique <- min(object$y):max(object$y)
        nUnique <- length(yUnique)
        p <- matrix(NA,length(yhat),nUnique)
        dimnames(p) <- list(NULL,yUnique)
        if(dist=="poisson"){
            p[,"0"] <- phi + (1-phi)*exp(-mu)
            for(i in 2:nUnique){
                p[,i] <- (1-phi)*dpois(x=yUnique[i],mu)
            }
        }
        if(dist=="negbin"){
            p[,"0"] <- phi + (1-phi)*((theta/(theta+mu))^theta)
            for(i in 2:nUnique){
                p[,i] <- (1-phi)*dnbinom(mu=mu,size=theta,x=yUnique[i])
            }
        }
        out$prob <- p
    }

    out$yhat <- yhat
    out$mu <- mu
    out$phi <- phi
   
    out
}

