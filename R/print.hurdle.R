print.hurdle <- function(x,digits=max(3,getOption("digits")-3), ...){
    if(!inherits(x,"hurdle"))
        stop("print.hurdle only implemented for objects of class hurdle")

    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")

    if(x$kz > 0){
        cat("\nHurdle Model ")
        if(x$link=="logit")
            cat("was fit with a logit link\n")
        if(x$link=="probit")
            cat("was fit with a probit link\n")
        print.default(format(coef(x)[1:x$kz]),
                             digits=digits,
                             print.gap=2,
                             quote=FALSE)
    }
    else cat("No coefficients for hurdle model\n")
    
    if(x$kx > 0){
        if(x$dist=="poisson"){
            cat("Count Model (Poisson)\n")
            cat("Coefficients:\n")
            print.default(format(coef(x)[(x$kz+1):(x$kz+x$kx)]),
                          digits = digits,
                          print.gap=2,
                          quote=FALSE)
        }
        
        if(x$dist=="negbin"){
            cat("Count Model (Negative Binomial)\n")
            cat("Coefficients:\n")
            print(x$coefficients[(x$kz+x$kx)+1,],
                  digits = digits,
                  ##signif.stars = signif.stars, 
                  ...)
            cat(paste("\nTheta =",
                      round(x$theta,digits),
                      "\n"))
        }
        
    }
    else cat("No coefficients for count component of  model\n")
    
    cat(paste("Number of Observations:",length(x$y),"\n"))
    cat(paste("Log-Likelihood:\t  ",format(round(x$llh,digits)),"\n"))

    invisible(NULL)
}
