"vuong" <-
function(m1,m2,digits=getOption("digits"))
{

    ## get predicted probabilities for both models
    m1class <- class(m1)[1]
    m2class <- class(m2)[1]

    m1y <- m1$y
    m2y <- m2$y
    m1n <- length(m1y)
    m2n <- length(m2y)
    if(m1n==0 | m2n==0)
        stop("could not extract dependent variables from models")

    if(m1n != m2n)
        stop(paste("models appear to have different numbers of observations\n",
                   "model 1 has ",m1n," observations\n",
                   "model 2 has ",m2n," observations\n",
                   sep="")
             )

    if(any(m1y != m2y))
        stop(paste("models appear to have different values on dependent variables\n"))

    whichCol <- match(m1y,min(m1y):max(m1y))

    m1p <- rep(NA,m1n)
    m2p <- rep(NA,m2n)
    p1 <- predprob(m1)
    p2 <- predprob(m2)
    for(i in 1:m1n){
        m1p[i] <- p1[i,whichCol[i]]
        m2p[i] <- p2[i,whichCol[i]]
    }

    m <- log(m1p/m2p)

    bad <- is.na(m) + is.nan(m) + (m==Inf) + (m==-Inf)
    if(any(bad)){
        cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
        cat("dropping these cases, but proceed with caution\n")
    }
    
    mbar <- mean(m[!bad])
    s <- sd(m[!bad])

    v <- sqrt(sum(!bad))*mbar/s

    cat(paste("Vuong Non-Nested Hypothesis Test-Statistic:",
              signif(v,digits),
              "\n"))
    cat("(test-statistic is asymptotically distributed N(0,1) under the\n")
    cat(" null that the models are indistinguishible)\n")

    cat("in this case:\n")
    if(v>0)
        cat(paste("model1 > model2, with p-value",signif(1-pnorm(v),digits),"\n"))
    else
        cat(paste("model2 > model1, with p-value",signif(pnorm(v),digits),"\n"))

    invisible(NULL)
}

