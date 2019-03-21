#' Checks whether the stard/end pos in pattern overlaps with the
#' start/end pos of each row in table by at least almin and
#' returns positions in table that overlap.
#'
#' @param pattern IR query pattern sequence positions.
#' @param table IR reference table positions to check for overlap with pattern.
#' @param aln.size Minimum allowable overlap considered valid.
#'
#' @return Positions in \code{table} that overlap.


isOverlapping <- function(pattern, table, aln.size = 7){

  opos <- IRanges::findOverlaps(makeIR(pattern), makeIR(table),
                       minoverlap = aln.size) %>% (S4Vectors::subjectHits)
  return(opos)
}
