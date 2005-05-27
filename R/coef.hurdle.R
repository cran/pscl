coef.hurdle <- function(object,...){
    if(!inherits(object,"hurdle"))
        stop("coef.hurdle only for objects of class hurdle")
    object$par
}
