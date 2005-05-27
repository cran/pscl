coef.zeroinfl <- function(object,...){
    if(!inherits(object,"zeroinfl"))
        stop("coef.zeroinfl only for objects of class zeroinfl")
    object$par
}
