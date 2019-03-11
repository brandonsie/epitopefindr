#' tidyBLAST
#'
#' Call a series of functions to tidy raw BLAST alignment table for use.
#'
#' @param blast Initial BLAST table to tidy.
#' @param fasta Path to fasta sequences to reference to annotate BLAST table
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#' with peptide sequences.

tidyBLAST <- function(blast, fasta, aln.size){
  blast %<>% organizeBLAST() #output table housekeeping (column names, etc.)
  blast %<>% numAlignments() #Add number of alignments per peptide
  blast %<>% addPepSeq(fasta) #add amino acid sequences (tile & align)
  blast %<>% decipherGaps() #split gapped alignments into smaller ungapped
  blast %<>% removeSmallAln(aln.size) #remove alignmens smaller than aln.size
  return(blast)
}
