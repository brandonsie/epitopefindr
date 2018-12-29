#' mergeFastaDuplicates
#' If two fasta entries have the same sequence, keep only one copy of the
#' sequence but concatenante the name to be Name1__Name2; __ can be replaced
#' with delimiter of choice.
#'
#' @param input XStringSet object to collapse.
#' @param sep Delimiter to place between names.

mergeFastaDuplicates <- function(input, sep = "__"){

  if(class(input)[1] == "AAStringSet"){
    input <- data.frame(ID = names(input), Seq = as.character(input))
  }

  input <- unique(input)
  if(sum(duplicated(input$Seq))>0){
    dup <- input$Seq[duplicated(input$Seq)] %>% unique
    for(i in 1:length(dup)){
      input$ID[input$Seq == dup[i]] %<>% paste(collapse=sep)}
    input <- input[!duplicated(input$Seq), ]
  }
  rownames(input) <- c(1:nrow(input))
  return(input)
}
