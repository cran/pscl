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
  
  ## gather up degrees of freedom
  k1 <- length(coef(m1))
  k2 <- length(coef(m2))
  
  m <- log(m1p) - log(m2p)  ## vector of log likelihood ratios (diffs of log probabilities)
  
  bad <- is.na(m) | is.nan(m) | is.infinite(m)
  neff <- sum(!bad)
  if(any(bad)){
    cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
    cat(paste("dropping these",sum(bad),"cases, but proceed with caution\n"))
  }
  
  aic.factor <- (k1-k2)/neff          
  bic.factor <- (k1-k2)/2 * log(neff)   
  print(bic.factor)
  
  ## test statistics
  v <- rep(NA,3)
  L1 <- sum(log(m1p[!bad]))
  L2 <- sum(log(m2p[!bad]))
  num <- rep(L1-L2,3) - c(0,aic.factor,bic.factor) 
  s <- sd(m[!bad])
  v <- num/(s*sqrt(neff))  
  names(v) <- c("Raw","AIC-corrected","BIC-corrected")
  print(v)
  print(s)
  print(num)
  
  ## bundle up for output
  pval <- rep(NA,3)
  msg <- rep("",3)
  for(j in 1:3){
    if(v[j]>0){
      pval[j] <- 1 - pnorm(v[j])
      msg[j] <- "model1 > model2"
    } else {
      pval[j] <- pnorm(v[j])
      msg[j] <- "model2 > model1"
    }
  }
  out <- data.frame(v,msg,format.pval(pval))
  names(out) <- c("Vuong z-statistic","H_A","p-value")

  ## output
  cat(paste("Vuong Non-Nested Hypothesis Test-Statistic:",
            "\n"))
  cat("(test-statistic is asymptotically distributed N(0,1) under the\n")
  cat(" null that the models are indistinguishible)\n")
  cat("-------------------------------------------------------------\n")
  
  print(out)

  return(invisible(NULL))
}

