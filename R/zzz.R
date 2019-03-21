#' .onload
#' Set stringsAsFactors FALSE for epitopefindr

#Initialize Envrionment
epitopefindrEnv <- new.env()

.onload <- function(libname, pkgname){
  options(stringsAsFactors = FALSE)

}