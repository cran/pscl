"predprob.zeroinfl" <-
function(obj,...){
    if(!inherits(obj,"zeroinfl"))
        stop("predprob.zeroinfl only available for zeroinfl objects\n")
    predict.zeroinfl(obj,type="prob",...)$prob
}

