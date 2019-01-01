#' qsSwap
#' Swap query and subject data from a BLAST table.
#' @param blast.in Input blast table to be swapped.
#' @export

qsSwap <- function(blast.in){
  #to force blast table to be symmetrical, swap subject and query info
  #make a list of which columns refer to query, subject, and either
  blast.in <- as.data.frame(blast.in)
  sNew <- qOld <- grep("^q|^Q",names(blast.in))
  qNew <- sOld <- grep("^s|^S",names(blast.in))

  blast.out <- data.frame(matrix(ncol=ncol(blast.in),nrow=nrow(blast.in))) %>%
    data.table::setnames(names(blast.in))

  for(i in 1:ncol(blast.in)){
    if(i %in% qNew){ #if this column should now be query
      qNewPos <- grep(paste0("^",i,"$"),qNew) #find which variable correspond
      qOldPos <- qOld[qNewPos] #and take that variable from the old position
      blast.out[,i] <- blast.in[,qOldPos]
    } else if(i %in% sNew){
      sNewPos <- grep(paste0("^",i,"$"),sNew)
      sOldPos <- sOld[sNewPos]
      blast.out[,i] <- blast.in[,sOldPos]
    } else{
      blast.out[,i] <- blast.in[,i]
    }
  }

  blast.out <- data.table::data.table(blast.out)

}
