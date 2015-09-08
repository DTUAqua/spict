#' @useDynLib spict
.onUnload <- function (lib) {
  library.dynam.unload("spict", lib)
}

.onAttach <- function(lib, pkg) {
   packageStartupMessage("Welcome to a Surplus Production model in Continuous Time (SPiCT).\n")  
 }
