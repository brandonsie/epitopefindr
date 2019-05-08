#' Call series of functions to standardize BLAST table for computation.
#'
#' @param blast Input tidied BLAST table.
#' @param fasta AAstringset object of fasta sequences.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#' @param tofilter Logical whether or not to call filterBLAST, which is only
#' necessary during initial setup.
#' @param use.doParallel Logical whether or not to use doParallel parallelization.
#'
#' @export

prepareBLAST <- function(blast, fasta, aln.size, use.doParallel = FALSE, tofilter=TRUE){
  #
  blast.merge <- rbind(blast, qsSwap(blast)) %>% unique
  if(nrow(blast.merge)==0){return(blast.merge)}
  blast.tidy <- tidyBLAST(blast.merge, fasta, aln.size, use.doParallel) #update col names, <7aa, gaps

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
