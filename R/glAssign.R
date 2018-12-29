#' glAssign
#' Quick wrapper function for assigning data to global environment
#'
#' @param name Identifier for object.
#' @param data Corresponding information.

glAssign <- function(name, data){

  #assign
  assign(paste0("gl.", name), data, envir = epitopefindrEnv)


}
