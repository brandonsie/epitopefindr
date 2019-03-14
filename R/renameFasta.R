#' renameFasta
#'
#' Revert peptide names of a fasta file back to original names based on naming map.
#'
#' @param fasta Blast table of alignments.
#' @param name.map Map of original and modified peptide names to convert.
#'
#' @export


renameFasta <- function(fasta, name.map){

  original <- name.map[,1]
  modified <- name.map[,2]

  fasta %<>% unmergeFastaDuplicates

  parsed.names <- names(fasta) %>% strsplit("\\.") %>% unlist %>%
    matrix(nrow = 3)

  basenames <- parsed.names %>% extract(1,)
  start <- parsed.names %>% extract(2,)
  end <- parsed.names %>% extract(3,)

  names(fasta) <- paste(original[match(basenames, modified)], start, end, sep = ".")

  return(fasta)
}