#' isOverlapping
#'
#' Checks whether the stard/end pos in pattern overlaps with the
#' start/end pos of each row in table by at least almin and
#' returns positions in table that overlap.
#'
#' @param pattern IR query pattern sequence positions.
#' @param table IR reference table positions to check for overlap with pattern.
#' @param almin Minimum allowable overlap considered valid.
#'
#' @return Positions in \code{table} that overlap.


isOverlapping <- function(pattern, table, almin = 7){

  opos <- IRanges::findOverlaps(makeIR(pattern), makeIR(table),
                       minoverlap = almin) %>% (S4Vectors::subjectHits)
  return(opos)
}
