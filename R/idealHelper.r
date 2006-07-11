## ideal helper functions

checkStart <- function(x, start) {
  if (as.numeric(start)>=(x$x[nrow(x$x),1]))
    stop("start must be less than the number of iterations")
  if(!as.integer(start)%in%x$x[,1])
    stop("start is not a valid iteration number")
  return (row(x$x)[x$x[,1]==as.integer(start),1])

}

checkD <- function(x,d) {
  if ((d<1)||(d>x$d))
    stop("d must be equal to one of the dimensions in the roll call object")
}

checkCI <- function(conf.int) {
  if((conf.int<=0)||(conf.int>=1))
      stop("conf.int must be between 0 and 1")
}

getDimX <- function(x,d,columns=TRUE) {
  checkD(x,d)
  px <- NULL
  if(columns) {
    px <- (x$x[,seq(from=d+1,to=ncol(x$x),by=x$d)])
    colnames(px)<-x$legis.names
  } else {
    px <- (x$x[seq(from=d,to=nrow(x$x),by=x$d),])
    rownames(px)<-x$legis.names
  }
  px
}

getDim <- function(x,d,dims,names,columns=TRUE) {
  px <- NULL
  if(columns) {
    px<-(x[,seq(from=d,to=ncol(x),by=dims)])
    colnames(px)<-names
  } else {
    px <- (x[seq(from=d,to=nrow(x),by=dims),])
    rownames(px)<-names
  }
  px
}
