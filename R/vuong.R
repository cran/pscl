"vuong" <- function(m1,m2,digits=getOption("digits")){
  ## get predicted probabilities for both models
  m1y <- m1$y
  m2y <- m2$y
  m1n <- length(m1y)
  m2n <- length(m2y)
  if(m1n==0 | m2n==0)
    stop("Could not extract dependent variables from models.")
  
  if(m1n != m2n)
    stop(paste("Models appear to have different numbers of observations.\n",
               "Model 1 has ",m1n," observations.\n",
               "Model 2 has ",m2n," observations.\n",
               sep="")
         )
  
  if(any(m1y != m2y)){
      stop(paste("Models appear to have different values on dependent variables.\n"))
  }

  whichCol <- match(m1y,min(m1y):max(m1y))  ## which column, matrix of predicted probs
  
  m1p <- rep(NA,m1n)
  m2p <- rep(NA,m2n)
  p1 <- predprob(m1)   ## likelihood contributions, model 1, cond on MLEs
  p2 <- predprob(m2)   ## likelihood contributions, model 2
  for(i in 1:m1n){
    m1p[i] <- p1[i,whichCol[i]]  ## pick off correct column, given observed y
    m2p[i] <- p2[i,whichCol[i]]
  }
  
  m <- log(m1p) - log(m2p)  ## vector of log likelihood ratios (diffs of log probabilities)
  
  bad <- is.na(m) | is.nan(m) | is.infinite(m)
  neff <- sum(!bad)
  if(any(bad)){
    cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
    cat("dropping these cases, but proceed with caution\n")
  }

  ## gather up degrees of freedom
  k1 <- length(coef(m1))
  k2 <- length(coef(m2))

  ## test statistic
  msum <- sum(m[!bad])
  s <- sd(m[!bad])
  v <- msum/(s * sqrt(neff))
  adj <- log(neff)*(k1/2 - k2/2)   ## adjustment a la AIC, length of model(s)
  v <- v - adj

  ## bundle up for output
  cat(paste("Vuong Non-Nested Hypothesis Test-Statistic:",
            signif(v,digits),
            "\n"))
  cat("(test-statistic is asymptotically distributed N(0,1) under the\n")
  cat(" null that the models are indistinguishible)\n")

  if(v>0){
      pval <- 1 - pnorm(v)
  } else {
      pval <- pnorm(v)
  } 
  
  cat("in this case:\n")
  if(v>0)
    cat(paste("model1 > model2, with p-value",
              format.pval(pval),
              "\n"))
  else
    cat(paste("model2 > model1, with p-value",
              format.pval(pval),
              "\n"))
  
  invisible(NULL)
}

