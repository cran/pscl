## parse and execute drop list used by summary.rollcall and ideal
dropRollCall <- function(object,dropList=NULL){
  if(class(object)!="rollcall"){
    stop("dropRollCall only works for objects of class rollcall.")
  }

  tmpRollCall <- object
  v <- tmpRollCall$votes
  
  if(!is.list(dropList)){
    cat("dropList must be a non-null list or alist. No subsetting will occur.")
    return(object)
  }

  cat("Dropping elements of rollcall matrix using the following dropList:\n")
  print(dropList)
  
  ## strip out user-designated votes of a particular code
  if(!is.null(dropList$codes) & length(dropList$codes)>0){
    cat("Processing dropList voting codes...\n")
    dc <- dropList$codes
    dCodes <- NULL
    if(all(is.character(dc))){   ## named element of codes list?
      dropCodes <- match(dc,names(tmpRollCall$codes))
      dropCodes <- dropCodes[!is.na(dropCodes)]
      #cat(paste("dropRollCall: dropCodes=",dropCodes,"\n"))

      if(length(dropCodes)>0){
        for(j in dropCodes)
          dCodes <- c(dCodes,tmpRollCall$codes[j])  ## drop these
        #cat(paste("dropRollCall: dCodes=",dCodes,"\n"))
        keepCodes <- !(names(tmpRollCall$codes) %in% dc)
        #cat(paste("dropRollCall: keepCodes=",keepCodes,"\n"))
        keepCodes <- tmpRollCall$codes[keepCodes]
        #cat(paste("dropRollCall: keepCodes=",keepCodes,"\n"))
        tmpRollCall$codes <- keepCodes    
        #cat(paste("dropRollCall: tmpRollCall$codes:\n"))
        #print(tmpRollCall$codes)
      }
    }
    if(is.numeric(dc)){    ## or numeric elements
      dCodes <- dc[dc %in% unique(as.vector(v))]
    }

    bad <- v %in% dCodes
    cat(paste("Will set",sum(bad),"voting decisions to NA.\n"))
    tmpRollCall$votes[bad] <- NA
    rm(bad)
  }

  dropLegis <- rep(FALSE,dim(v)[1])
  dropVotes <- rep(FALSE,dim(v)[2])

  ## drop legislators if too little data
  if(!is.null(dropList$legisMin)){
    legisMin <- dropList$legisMin
    if(length(legisMin)!=1 |
       is.na(legisMin) |
       !is.numeric(legisMin) |
       legisMin >= tmpRollCall$m)
      stop("bad value for legisMin in drop list.")
    vtmp <- convertCodes(tmpRollCall)
    goodCount <- apply(vtmp,1,function(x)sum(!is.na(x)))
    dropLegis <- dropLegis | goodCount<legisMin
  }

  ## check for subsetting in legis.data
  if(!is.null(dropList$dropLegis)){
    r <- dropRollCallViaData(dropList$dropLegis,
                             object=tmpRollCall,
                             d=expression(legis.data))
    if(!is.null(r))
      dropLegis <- dropLegis | r
  }
  
  ## drop votes by lop-sidedness
  if(!is.null(dropList$lop)){
    cat("Processing lop-sided restrictions...\n")
    lop <- dropList$lop
    if(length(lop)!=1 |
       is.na(lop) |
       !is.numeric(lop) |
       lop < 0 |
       lop >= tmpRollCall$n)
      stop("Invalid value for lop")
    if(is.null(tmpRollCall$voteMargins)){
      cat("Computing vote margins...\n")
      tmpRollCall <- computeMargins(tmpRollCall,dropList=NULL)
    }
    r <- tmpRollCall$voteMargins[,"Min"] <= lop
    cat(paste("Will drop",sum(r),"roll calls.\n"))
    dropVotes <- dropVotes | r
    cat("Finished processing lop-sided restrictions.\n")
  }

  ## check for subsetting in vote.data
  if(!is.null(dropList$dropVotes)){
    r <- dropRollCallViaData(dropList$dropVotes,
                             object=tmpRollCall,
                             d=expression(vote.data))
    if(!is.null(r))
      dropVotes <- dropVotes | r
  }

  ## final processing
  cat(paste("dropRollCall will drop",
            sum(dropLegis),
            "legislators and",
            sum(dropVotes),
            "rollcalls.\n"))
  if(sum(dropLegis)>0){
    cat("Dropped Legislators:\n")
    print(dimnames(tmpRollCall$votes)[[1]][dropLegis])
  }
  if(sum(dropVotes)>0){
    cat("Dropped Votes:\n")
    print(dimnames(tmpRollCall$votes)[[2]][dropVotes])
  }
    
  tmpRollCall$votes <- tmpRollCall$votes[!dropLegis,!dropVotes]

  if(!is.null(tmpRollCall$legis.data))
    tmpRollCall$legis.data <- tmpRollCall$legis.data[!dropLegis,]

  if(!is.null(tmpRollCall$vote.data))
    tmpRollCall$vote.data <- tmpRollCall$vote.data[!dropVotes,]

  if(!is.null(tmpRollCall$voteMargins))
    tmpRollCall$voteMargins <- tmpRollCall$voteMargins[!dropVotes,]

  tmpRollCall$n <- dim(tmpRollCall$votes)[1]
  tmpRollCall$m <- dim(tmpRollCall$votes)[2]

  tmpRollCall                      ## return rollcall object
}

dropRollCallViaData <- function(expr,object,d){
  cf <- match.call()

  f <- try(eval(d,envir=object),silent=TRUE)
  if(inherits(f,"try-error")){
    cat(paste("The data frame ",
              cf$d,
              " was not found in ",
              cf$object,
              ".\n",sep=""))
    cat("Proceeding by ignoring this subsetting restriction.\n")
    return(NULL)
  }

  r <- try(eval(expr,f),silent=TRUE)
  if(inherits(r,"try-error")){
    r <- rep(FALSE,dim(f)[1])
    cat(paste("The assertion ",
              deparse(expr),
              " could not be evaluated in the ",
              cf$d,
              " component of ",
              cf$object,
              ".\n",
              sep=""))
    cat("Proceeding by ignoring this assertion.\n")
  }
  if(!is.logical(r))
    stop("'x' must evaluate to logical")
  r <- r & !is.na(r)
  r
}
