#' addPepSeq
#' Add/update column in blast table listing peptide names.
#'
#' @param blast BLAST table to update with sequences.
#' @param fasta Fasta file with relevant peptide names and amino acid sequences.

addPepSeq <- function(blast, fasta){
  if(class(fasta) == "character") fasta <- Biostrings::readAAStringSet(fasta)
  fdf <- data.frame(id = names(fasta), seq = fasta %>% as.character,
                    stringsAsFactors = FALSE)

  #add/update column in blast table listing peptide sequence
  for(i in 1:nrow(blast)){
    blast$qSeq[i] <- fdf$seq[fdf$id == blast$qID[i]]
    blast$sSeq[i] <- fdf$seq[fdf$id	 == blast$sID[i]]
  }
  return(blast)
}
