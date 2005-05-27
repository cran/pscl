"logLik.hurdle" <-
    function(object,...){
        if(!inherits(object,"hurdle"))
            stop("this function only works for objects of class hurdle\n")
        object$llh
    }
