#' Remove filler sequences after "." stop codon and cterminal "*"
#'
#' @param input Untidy fasta.
#'
#' @export

tidyFasta <- function(input){

  if(class(input) == "character"){
    name <- input
    input <- Biostrings::readAAStringSet(input)
  }

  fasta.df <- data.frame(ID = (input %>% names %>% as.character),
                         Seq = (input %>% as.character))

  if(length(unique(fasta.df$ID)) < length(fasta.df$ID)){
    stop("Error: all input peptides should have unique names.")
  }

  if(grepl("__________", fasta.df$ID) %>% mean > 0){
    stop("Error: input peptide names must not contain 10 consecutive underscores.")
  }

  #convert non-alphanumeric characters to hyphens
  #then remove any leading/lagging hyphens and multi-hyphens
  fasta.df$ID %<>% gsub("[^[:alnum:]]", "-", .)

  while(length(grep("--", fasta.df$ID)) > 0){
    fasta.df$ID %<>% gsub("--","-", .)
  }
  fasta.df$ID %<>% gsub("^-", "", .)
  fasta.df$ID %<>% gsub("-$", "", .)


  fasta.df$ID %<>% gsub("^-", "", .)



  if(length(unique(fasta.df$ID)) < length(fasta.df$ID)){
    stop("Error: epitopefindr (temporarily) coerces non-alphanumeric characters in peptide names to hypens '-'. All peptides must have unique names when only alphanumeric characters are retained.")
  }



  fasta.df$Seq %<>% strsplit("\\.")
  fasta.df$Seq %<>% sapply(function(x){x[[1]][1]})
  fasta.df$Seq %<>% gsub("\\*", "", .)

  output <- Biostrings::AAStringSet(fasta.df$Seq)
  names(output) <- fasta.df$ID
  return(output)
}
