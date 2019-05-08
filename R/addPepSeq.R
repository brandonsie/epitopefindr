#' Add/update column in blast table listing peptide names.
#'
#' @param blast BLAST table to update with sequences.
#' @param fasta Fasta file with relevant peptide names and amino acid sequences.

addPepSeq <- function(blast, fasta){
  if(class(fasta) == "character"){fasta <- Biostrings::readAAStringSet(fasta)}
  fdf <- data.frame(id = names(fasta), seq = fasta %>% as.character,
                    stringsAsFactors = FALSE)

  #add/update column in blast table listing peptide sequence

  unique_peptides <- c(blast$qID, blast$sID) %>% unique
  blast$qSeq <- blast$sSeq <- NA

  for(i in 1:length(unique_peptides)){
    q_rows <- c(1:nrow(blast))[blast$qID == unique_peptides[i]]
    if(length(q_rows) > 0) {blast$qSeq[q_rows] <-
      fdf$seq[fdf$id == unique_peptides[i]]}

    s_rows <- c(1:nrow(blast))[blast$sID == unique_peptides[i]]
    if(length(s_rows) > 0) {blast$sSeq[s_rows] <-
      fdf$seq[fdf$id == unique_peptides[i]]}
  }


  return(blast)
}
