#' Update BLAST table positions to reflect smaller intervals when a subinterval
#' has been determined as the minimal overlap of alignments to a peptide.
#'
#' @param data List containing BLAST table and fasta file and index peptide order to process.
#' @param tofilter Binary whether or not to filter BLAST table entries.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#'
#' @export

trimEpitopes <- function(data, aln.size, tofilter = FALSE){

  blast <- data[[1]]
  fasta <- data[[2]]
  index.order <- data[[3]]

  blast %<>% rbind(qsSwap(blast)) %>% unique
  blast %<>% removeSmallAln(aln.size)

  #update blast table in reverse order
  pb <- epPB(1,length(index.order))
  for(i in length(index.order):1){
    utils::setTxtProgressBar(pb,length(index.order)-(i-1))
    index <- index.order[i] %>% as.character

    # print(i)
    # print(index)
    blast.backup <- blast

    # Getting error in indexEpitope which seems to stem from double qsSwap() calls. qsSwap() has already been called above on line 16. The double call let too short queries that were previously subjects enter the algorithm, as removeSmallAln() on line 17 only operates on queries.
    blast.index <- 
      blast[blast$qID == index, -"nAlign"] %>% 
      unique %>% 
      numAlignments()
    

    if(nrow(blast.index)>0){
      #input full blast and index name. output modified blast and index epitopes
      indexData <- indexEpitopes(blast, index, aln.size)
      blast <- indexData[[1]]
    }
  }
  close(pb)


  #output
  # mpath <- "final_epitopes.fasta"
  finalep <- blast[order(blast$qID,blast$qStart,blast$qEnd),
                        c("qID", "qStart", "qEnd", "qSeq")] %>% unique
  finalep$Seq <- sapply(1:nrow(finalep), function(x){
    substr(finalep$qSeq[x], finalep$qStart[x], finalep$qEnd[x])
  })
  finalep$ID <- paste(finalep$qID, finalep$qStart, finalep$qEnd, sep=".")
  finalep %<>% mergeFastaDuplicates
  print(paste(nrow(finalep),"epitope sequences identified."))
  # writeFastaAA(finalep, mpath)


  final.stringset <- Biostrings::AAStringSet(finalep$Seq)
  names(final.stringset) <- finalep$ID

  # writeFastaAA(final.stringset,"testfinalep.fasta")

  outputData <- list(blast = blast, fasta = final.stringset)

  return(outputData)
}
