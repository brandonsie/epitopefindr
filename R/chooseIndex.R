#` chooseIndex
#' Identifies which of a set of remaining peptides from a blast
#' alignment table should serve as the next index peptide.
#'
#' @param input Data frame with BLAST alignments.
#' @param method Either "min" or "max": instructs whether chooseIndex should
#' suggest the peptide with the fewest or most alignments, respectively.
#' @return Peptide name.

chooseIndex <- function(input, method="min"){


  #list of remaining full peptide sequences
  fullpep <- input[!grepl("\\.", input$qID), c("qID","nAlign")] %>% unique

  if(method == "min"){ #by default, select index peptide with FEWEST alignments
    index <- fullpep$qID[fullpep$nAlign == min(fullpep$nAlign)][1] %>%
      as.character
  } else if(method == "max"){
    index <- fullpep$qID[fullpep$nAlign == max(fullpep$nAlign)][1] %>%
      as.character
  } else{stop("Error: chooseIndex: improper selection method.")}
}

