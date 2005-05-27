odTest <- 
function(glmobj,digits=max(3,getOption("digits")-3)){
    if(class(glmobj)[1]!="negbin")
        stop("this function only works for objects of class negbin\n")
    
    poissonGLM <- glm(formula=eval(glmobj$call$formula),
                      data=eval(glmobj$call$data),
                      family="poisson")
    require(stats)
    llhPoisson <- logLik(poissonGLM)
    llhNB <- logLik(glmobj)
    
    d <- 2*(llhNB - llhPoisson)
    pval <- 1 - pchisq(d,df=1)
    
    cat("Likelihood ratio test of HO: Poisson, as restricted NB model:\n")
    cat(paste("Chi-Square Test Statistic = ",
              round(d,digits),
              "p-value =",
              format.pval(pval,digits=digits),
              "\n"))

    invisible(NULL)
}
