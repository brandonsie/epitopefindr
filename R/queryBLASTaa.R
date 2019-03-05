#' queryBLASTaa
#'
#' BLAST query sequence against a database of sequences.
#'
#' @param query Path to fasta file of query sequence.
#' @param db Path to fasta file of database sequences.
#'
#' @export

queryBLASTaa <- function(query, db){
  #
  rBLAST::makeblastdb(db, dbtype="prot")
  blastdb <- rBLAST::blast(db, type="blastp")
  blastseq <- Biostrings::readAAStringSet(query)
  blastpred <- stats::predict(blastdb, blastseq)
}
