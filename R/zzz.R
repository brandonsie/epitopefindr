#' .onload
#' Set stringsAsFactors FALSE for epitopefindr

#Initialize Envrionment
epitopefindrEnv <- new.env()

.onload <- function(libname, pkgname){
  options(stringsAsFactors = FALSE)





  # data.frame <- function (..., stringsAsFactors = FALSE)
  # {
  #   base::data.frame(..., stringsAsFactors = stringsAsFactors)
  # }
  #
  # as.data.frame <- function (x, row.names = NULL, optional = FALSE, stringsAsFactors = FALSE,
  #           ...)
  # {
  #   base::as.data.frame(x, row.names = NULL, optional = FALSE,
  #                       stringsAsFactors = stringsAsFactors, ...)
  # }

}