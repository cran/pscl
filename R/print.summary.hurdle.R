"print.summary.hurdle" <-
function(x,digits = max(3, getOption("digits") - 3),...)
{                               
  cat("Hurdle Model Summary\n")
  cat("--------------------\n")

  cat("\nCall:\n")
  cat(paste(deparse(x$call),sep="\n",collapse="\n"),
      "\n\n",sep="")

  cat(paste("Total Log-likelihood:",x$llh,"\n"))
  
  cat("\nHurdle Model (Zeros vs Other Outcomes) ")
  if(x$link=="logit")
    cat("was fit by logit")
  if(x$link=="probit")
    cat("was fit by probit")

  cat("\nCoefficients:\n") 
  print(x$coefficients[1:x$kz,],
        digits = digits, 
        ...)

  if(x$dist=="poisson"){
    cat("\n\nCount Model for Non-Zero Outcomes (Poisson)\n")
    cat("Coefficients:\n")
    print(x$coefficients[(x$kz+1):(x$kz+x$kx),],
          digits = digits,
          ...)
  }

  if(x$dist=="negbin"){
    cat("\n\nCount Model for Non-Zero Outcomes (Negative Binomial)\n")
    cat("Coefficients:\n")
    print(x$coefficients[(x$kz+1):(x$kz+x$kx+1),],
          digits = digits,
          ...)
  }
  invisible(NULL)
}

