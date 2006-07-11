"summary.hurdle" <-
function(object, ...){
  out <- list()
  out$call <- object$call
  out$y <- object$y
  out$n <- object$n
  out$kz <- object$kz
  out$kx <- object$kx
  out$dist <- object$dist
  out$gamma <- object$par[1:object$kz]
  out$beta <- object$par[(object$kz+1):(object$kz+object$kx)]
  if(object$dist=="negbin")
    out$theta <- object$theta
  out$vc <- -solve(object$hessian)
  se <- sqrt(diag(out$vc))
  out$coefficients <- cbind(object$par,
                            se)
  tstat <- out$coefficients[,1]/out$coefficients[,2]
  pval <- 2*pnorm(-abs(tstat))
  out$coefficients <- cbind(out$coefficients,tstat,pval)
  dimnames(out$coefficients) <- list(names(object$par),
                                     c("Estimate",
                                       "Std. Error",
                                       "z value",
                                       "Pr(>|z|)"))
  out$llh <- object$llh
  out$method <- object$method
  out$link <- object$link
  out$call <- object$call
  
  class(out) <- "summary.hurdle"
  out
}

