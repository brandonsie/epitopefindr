#' removeSmallAln
#'
#' Delete rows in blast table for which the alignment is fewer than 7aa or
#' other custom threshold.
#'
#' @param blast Input BLAST table.
#' @param minsize Minimum length of alignment to keep.
#'
#' @export

removeSmallAln <- function(blast, minsize = 7){
  #
  toosmall <- (1:nrow(blast))[(blast$qEnd-blast$qStart) < minsize]
  if(length(toosmall)>0){blast <- blast[-toosmall, ]}
  return(blast)
}
