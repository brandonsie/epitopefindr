#' glGet
#'
#' Retrieve data from global environment, e.g. from previous glAssign call.
#'
#' @param name Varable to retrieve.

glGet <- function(name){

  get(paste0("gl.", name), envir = epitopefindrEnv)
}
