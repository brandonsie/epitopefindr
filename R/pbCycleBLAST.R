#' pbCycleBLAST
#' Initialize a txtProgressBar object with which to call cycleBLAST.
#' @param path Path to current fasta file to pass to cycleBLAST.
#' @param ncycles Number of cycles of cycleBLAST to perform.
#' @export

pbCycleBLAST <- function(path, ncycles="max"){
  # == == == == == A. Initialize == == == == ==
  blast.main <- glGet("blast.main")
  path <- glParamGet("path")

  # == == == == == B. Count starting full peptides. == == == == ==
  if(ncycles == "max"){
    n.pep <- blast.main$qID %>% unique %>% length

  } else if(class(ncycles) == "numeric"){
    n.pep <- ncycles

  } else {
    stop("Error: pbCycleBLAST: improper ncycles parameter.")
  }

  # == == == == == C. Run cycle of epitopeBLAST. == == == == ==
  pb <- epPB(-n.pep,0)

  path %<>% cycleBLAST(pb, n.pep)
  glParamAssign("path", path)

  return(path)

}
