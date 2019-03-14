#' renameBLAST
#'
#' Revert qID and sID columns of blast table back to original names based on naming map.
#'
#' @param blast Blast table of alignments.
#' @param name.map Map of original and modified peptide names to convert.
#'
#' @export


renameBLAST <- function(blast, name.map){

  original <- name.map[,1]
  modified <- name.map[,2]

  blast$qID <- original[match(blast$qID, modified)]
  blast$sID <- original[match(blast$sID, modified)]

  return(blast)
}