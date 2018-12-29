#' trimEpitopes
#' Update BLAST table positions to reflect smaller intervals when a subinterval
#' has been determined as the minimal overlap of alignments to a peptide.
#' @param path Path to current fasta file.
#' @param tofilter Binary whether or not to filter BLAST table entries.
#' @export

trimEpitopes <- function(path, tofilter = FALSE){

  #load some data from global environment
  gl <- c("blast.id3","output.dir")
  for(i in gl){assign(i,glParamGet(i))}

  blast.main <- data.table::fread(blast.id3)
  blast.main %<>% prepareBLAST(tofilter)
 glAssign("blast.main", blast.main)

  index.order <- fread(paste0(output.dir,"epitopes/indexOrder.txt"),
                       header = FALSE)



  #update blast table in reverse order
  pb <- epPB(1,nrow(index.order))
  for(i in nrow(index.order):1){
    utils::setTxtProgressBar(pb,nrow(index.order)-i)
    blast.main <- data.table::as.data.table(glGet("blast.main"))
    index <- index.order[i] %>% as.character

    # print(i)
    # print(index)

    blast.index <- rbind(
      blast.main[blast.main$qID==index,-"nAlign"],
      qsSwap(blast.main[blast.main$sID==index,-"nAlign"])) %>% unique
    blast.backup <- blast.main

    if(nrow(blast.index)>0){
      indexEpitopes(as.data.frame(blast.index))
    }
  }
  close(pb)


  blast.main <- glGet("blast.main")

  #output
  mpath <- paste0(output.dir, "epitopes/final_epitopes.fasta")
  finalep <- blast.main[order(blast.main$qID,blast.main$qStart,blast.main$qEnd),
                        c("qID", "qStart", "qEnd", "qSeq")] %>% unique
  finalep$Seq <- sapply(1:nrow(finalep), function(x){
    substr(finalep$qSeq[x], finalep$qStart[x], finalep$qEnd[x])
  })
  finalep$ID <- paste(finalep$qID, finalep$qStart, finalep$qEnd, sep=".")
  writeFastaAA(finalep %>% mergeFastaDuplicates, mpath)
  print(paste(nrow(finalep),"epitope sequences identified."))

  glParamAssign("path", mpath)
  return(mpath)
}
