#' mergeFastaDuplicates
#' If two fasta entries have the same sequence, keep only one copy of the
#' sequence but concatenante the name to be Name1__Name2; __ can be replaced
#' with delimiter of choice.
#'
#' @param input XStringSet object to collapse.
#' @param sep Delimiter to place between names.
#' @export

mergeFastaDuplicates <- function(input, sep = "__"){

  options(stringsAsFactors = FALSE)
  if(class(input)[1] == "AAStringSet"){
    input <- data.frame(ID = names(input), Seq = as.character(input))
  }
  output <- unique(input)
  output$ID %<>% as.character
  output$Seq %<>% as.character

  if(sum(duplicated(output$Seq))>0){
    dup <- output$Seq[duplicated(output$Seq)] %>% unique
    for(i in 1:length(dup)){
      output$ID[output$Seq == dup[i]] %<>% paste(collapse=sep)}
    output <- output[!duplicated(output$Seq), ]
  }
  rownames(output) <- c(1:nrow(output))
  return(output)
}
