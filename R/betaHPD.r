betaHPD <- function(alpha,beta,p=.95,plot=FALSE){

    if(is.na(p) | is.nan(p) | p > 1 | p < 0)
        stop("p not between 0 and 1\n")

    if(alpha<=1 | beta <=1)
        stop("betaHPD only implemented for alpha and beta both > 1\n")

    func <- function(x0,alpha,beta){
        y0 <- dbeta(x0,alpha,beta)
        p0 <- pbeta(x0,alpha,beta)
        x1 <- qbeta(p0+p,alpha,beta)
        y1 <- dbeta(x1,alpha,beta)
        out <- abs(y0-y1)
        out
    }
    
    foo <- try(optimize(f=func,alpha=alpha,beta=beta,
                        interval=c(.Machine$double.eps,
                        qbeta(1-p,
                              alpha,beta))))
    if(inherits(foo,"try-error")){
        warning("optimization in betaHPD failed\n")
        out <- rep(NA,2)
    }
    else{
        out <- c(foo$minimum,
                 qbeta(pbeta(foo$minimum,alpha,beta)+p,
                       alpha,beta)
                 )
        if(plot){
            xseq <- seq(min(qbeta(.0001,alpha,beta),out[1]),
                        max(qbeta(.9999,alpha,beta),out[2]),
                        length=1000)
            plot(xseq,dbeta(xseq,alpha,beta),
                 xlab=expression(theta),
                 ylab="",
                 axes=F,
                 type="n")
            axis(1)
            dseq <- seq(out[1],out[2],length=250)
            fx <- dbeta(dseq,alpha,beta)
            polygon(x=c(out[1],dseq,rev(dseq)),
                    y=c(0,fx,rep(0,250)),
                    border=F,col=gray(.45))
            lines(xseq,dbeta(xseq,alpha,beta))
        }
    }
    out
}
