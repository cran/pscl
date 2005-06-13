.onAttach <- function(...){
    cat(paste("   pscl",
              paste(rep(".",floor(getOption("width")*.90 - 4)),collapse=""),
              "\n",
              sep="")
        )
    cat("   R classes and methods developed in the\n")
    cat("   Political Science Computational Laboratory\n")
    cat("   Department of Political Science, Stanford University\n")
    cat("   Simon Jackman <jackman@stanford.edu>\n")
    cat("   http://pscl.stanford.edu\n")
    invisible(NULL)
}

.onUnload <- function(){
    invisible(NULL)
}
