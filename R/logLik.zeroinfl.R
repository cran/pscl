"logLik.zeroinfl" <-
    function(object,...){
        if(!inherits(object,"zeroinfl"))
            stop("this function only works for objects of class zeroinfl\n")
        object$llh
    }
