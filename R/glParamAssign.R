#' glParamAssign
#' Special case of glAssign to update a parameter in a table called gl.params.
#'
#' @param name ID of variable to populate first column "param" of gl.params
#' @param data ID of value to populate second column "value" of gl.params

glParamAssign <- function(name,data){

  #assign
  if(!("gl.params" %in% ls(epitopefindrEnv))){
    assign("gl.params", matrix(nrow=0,ncol=2) %>% data.frame %>%
             data.table::setnames(c("param","value")), envir = epitopefindrEnv)
  }

  params <- get("gl.params", envir = epitopefindrEnv)

  #update value if exists. otherwise create new
  if(name %in% params$param){
    params$value[params$param == name] <- data
  } else{
    params <- rbind(params,data.table::setnames(data.frame(t(c(name,data))),names(params)))
  }

  assign("gl.params", params, envir = epitopefindrEnv)



}
