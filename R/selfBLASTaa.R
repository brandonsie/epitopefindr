#' selfBLASTaa
#' Run BLASTp on amino acid sequences against themselves.
#'
#' @param path Path to fasta file of amino acid sequences to BLAST.
#' @export

selfBLASTaa <- function(path){
  rBLAST::makeblastdb(path, dbtype="prot")
  blastdb <- rBLAST::blast(path, type="blastp")
  blastseq <- Biostrings::readAAStringSet(path)
  blastpred <- stats::predict(blastdb, blastseq)
}
