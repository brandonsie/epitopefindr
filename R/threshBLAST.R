#' threshBLAST
#' Subset initial BLAST alignbment table to only include alignments below a
#' specified e-value threshold.
#' @param blast Initial BLAST alignment table.
#' @param e.thresh Maximum e-value threshold to retain.
#' @export

threshBLAST <- function(blast, e.thresh){
  blast <- blast[blast$E < as.numeric(e.thresh), ]
  if(nrow(blast) > 0){return(blast)}
}