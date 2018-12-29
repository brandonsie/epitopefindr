#'tidyFasta
#' remove filler sequences after "." stop codon and cterminal "*"
#'
#' @param input Untidy fasta.
#' @param name Filename to write tidy fasta file.

tidyFasta <- function(input, name){

  if(class(input) == "character"){
    name <- input
    input <- Biostrings::readAAStringSet(input)
  }

  fasta.df <- data.frame(ID = (input %>% names %>% as.character),
                         Seq = (input %>% as.character))

  fasta.df$ID %<>% gsub(" ", "-", .)
  fasta.df$ID %<>% gsub("\\.", "-", .)
  fasta.df$ID %<>% gsub("_", "-", .)
  fasta.df$ID %<>% gsub("\\(|\\)","",.)
  fasta.df$ID %<>% gsub(",","",.)

  fasta.df$Seq %<>% strsplit("\\.")
  fasta.df$Seq %<>% sapply(function(x){x[[1]][1]})

  fasta.df$Seq %<>% gsub("\\*", "", .)

  seqinr::write.fasta(as.list(fasta.df$Seq), fasta.df$ID, name)
  output <- Biostrings::readAAStringSet(name)
}
