#' numAlignments
#'
#' Add/update column to blast table listing # of alignments per query sequence.
#'
#' @param input BLAST table to update with number of alignments per peptide.

numAlignments <- function(input){
  input$nAlign <- apply(input, 1, function(x) sum(input$qID == x[1]))
  return(input)
}
