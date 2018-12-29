#' freadSF
#' Custom data.table::fread wrapper with stringsasfactors set to FALSE
#' @param input File to read into R.


freadSF <- function(input){
  data.table::fread(input,stringsAsFactors = FALSE)
}
