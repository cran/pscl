## ## for creating a package
.package.Name <- "pscl"

##.First.lib <- function(lib,pkg){
##  library.dynam(.package.Name,
##                pkg,
##                lib)
##}

##.Last.lib <- function(libpath){
##  library.dynam.unload(chname="pscl",libpath=libpath)
##}

.onAttach <- function(...){
  cat("Classes and Methods for R developed in the\n")
  cat("Political Science Computational Laboratory\n")
  cat("Department of Political Science\n")
  cat("Stanford University\n")
  cat("Simon Jackman\n")
  cat("hurdle and zeroinfl functions by Achim Zeileis\n") 
}

.onUnload <- function(libpath){
  library.dynam.unload("pscl",libpath=libpath)
}
