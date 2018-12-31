#' writeFastaAA
#' Takes data frame with $Seq and $ID or AAStringSet and writes fasta file
#'
#' @param input Data frame or AAStringSet to write to fasta file.
#' @param fpath Directory path to write fasta file to.
#' @export

writeFastaAA <- function(input, fpath){

  if(class(input)[1] == "data.table"){input <- data.frame(input)}

  if(class(input) == "data.frame"){
    seqinr::write.fasta(input$Seq %>% as.list, input$ID, fpath, nbchar = 90, as.string=TRUE)
  } else if(class(input) == "AAStringSet"){
    seqinr::write.fasta(input %>% as.character %>% as.vector %>% as.list,
                names(input), fpath, as.string=TRUE, nbchar = 90)
  } else stop("ERROR: unrecognized input to writeFastaAA")
}
