#' pbCycleBLAST
#' Initialize a txtProgressBar object with which to call cycleBLAST.
#' @param blast BLAST alignment table to process.
#' @param fasta Fasta file of peptide sequences to process.
#' @param ncycles Number of cycles of cycleBLAST to perform.
#' @export

pbCycleBLAST <- function(blast, fasta, ncycles="max"){
  # == == == == == A. Initialize == == == == ==
  options(stringsAsFactors = FALSE)
  # == == == == == B. Count starting full peptides. == == == == ==
  if(ncycles == "max"){
    n.pep <- blast$qID %>% unique %>% length
  } else if(class(ncycles) == "numeric"){
    n.pep <- ncycles
  } else {
    stop("Error: pbCycleBLAST: improper ncycles parameter.")
  }

  # == == == == == C. Run cycle of epitopeBLAST. == == == == ==
  pb <- epPB(-n.pep,0)

  data <- list(blast = blast, fasta = fasta)
  data %<>% cycleBLAST(pb, n.pep)

  return(data)

}
