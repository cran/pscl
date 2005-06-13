"print.summary.zeroinfl" <-
function(x,digits = max(3, getOption("digits") - 3),...)
{                               

    titleString <- "Zero-Inflated Count Model Summary"
    zeroModelString <- "Zero-Inflated Model"
    cat(titleString,"\n")
    cat("---------------------------\n")
    
    cat("\nCall:\n")
    cat(paste(deparse(x$call),sep="\n",collapse="\n"),
        "\n\n",sep="")

    cat(paste("Total Log-likelihood:",x$llh,"\n"))
    
    cat("\n",zeroModelString)
    if(x$link=="logit")
        cat(" was fit with a logit link\n")
    if(x$link=="probit")
        cat(" was fit with a probit link\n")

  kx <- x$kx
  kz <- x$kz
  
  cat("Coefficients:\n") 
  print.matrix(x$coefficients[1:kz,],
               digits = digits,
               ##signif.stars = signif.stars, 
                ...)
  cat("---------------------------------------------------------------------------\n")
  
  if(x$dist=="poisson"){
    cat("\nCount Model (Poisson)\n")
    cat("Coefficients:\n")
    print.matrix(x$coefficients[(kz+1):(kz+kx),],
                 digits = digits,
                 ##signif.stars = signif.stars, 
                  ...)
  }

  if(x$dist=="negbin"){
    cat("\nCount Model (Negative Binomial)\n")
    cat("Coefficients:\n")
    print.matrix(x$coefficients[(kz+1):((kz+kx)+1),],
                  digits = digits,
                  ##signif.stars = signif.stars, 
                  ...)
    cat(paste("\nTheta =",
              round(x$theta,digits),
              "\n"))
  }
  invisible(NULL)
}

