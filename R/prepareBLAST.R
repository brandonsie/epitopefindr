#' prepareBLAST
#' Call series of functions to standardize BLAST table for computation.
#' @param blast.thresh Input tidied BLAST table.
#' @param tofilter Logical whether or not to call blast.filter, which is only
#' necessary during initial setup.

prepareBLAST <- function(blast.thresh, tofilter=TRUE){
  #
  path0 <- glParamGet("path0")
  fasta <- Biostrings::readAAStringSet(path0)
  blast.merge <- rbind(blast.thresh, qsSwap(blast.thresh)) %>% unique
  if(nrow(blast.merge)==0){return(blast.merge)}
  blast.tidy <- tidyBLAST(blast.merge, fasta) #update col names, <7aa, gaps

  if(!exists("tofilter")){tofilter <- TRUE}
  if(tofilter){
    blast.filter <- filterBLAST(blast.tidy) #remove self-alignments
    if(nrow(blast.filter)==0){return(blast.filter)}
    rownames(blast.filter) <- c(1:nrow(blast.filter)) #fix row numbering
    return(blast.filter)
  } else{
    if(nrow(blast.tidy)==0){return(blast.tidy)}
    rownames(blast.tidy) <- c(1:nrow(blast.tidy))
    return(blast.tidy)
  }
}
