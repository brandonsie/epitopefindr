#' Delete rows in blast table for which the alignment is fewer than 7aa or
#' other custom threshold.
#'
#' @param blast Input BLAST table.
#' @param aln.size Minimum length of alignment to keep.
#'
#' @export

removeSmallAln <- function(blast, aln.size){
  #
  toosmall <- (1:nrow(blast))[(blast$qEnd-blast$qStart) < (aln.size - 1)]
  if(length(toosmall)>0){blast <- blast[-toosmall, ]}
  return(blast)
}
