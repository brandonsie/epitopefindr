#' epSetupPeptides
#'
#' Tidy input peptide sequences.

epSetupPeptides <- function(){
  gl <- c("path","output.dir","proj.id")
  for(i in gl){assign(i,glParamGet(i))}

  print(path)
  fasta <- tidyFasta(path) #remove filler .&* Update path & return AAString
  glAssign("fasta", fasta)

  fpath <- paste0(output.dir, "epitopes/", proj.id, ".fasta")
  writeFastaAA(fasta, fpath)
  glParamAssign("pathi",path)
  glParamAssign("path", fpath)
  glParamAssign("path0", fpath)
}
