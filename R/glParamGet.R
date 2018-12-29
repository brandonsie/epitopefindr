#' glParamGet
#' Special case of glGet to retrieve a value from object gl.params.
#'
#' @param name Variable to retrieve.

glParamGet <- function(name){
  params <- get("gl.params", envir = epitopefindrEnv)
  return(params$value[params$param == name])
}
