#' Split gapped BLAST alignments into smaller ungapped alignments.
#'
#' @param blast Input gapped BLAST alignment table to ungap.
#' @param use.doParallel Logical whether or not to use doParallel parallelization.

decipherGaps <- function(blast, use.doParallel){
  gpos <- c(1:nrow(blast))[blast$Gaps>0]
  if(length(gpos) == 0){return(blast)}


  if(!use.doParallel){
    pb <- epPB(0,length(gpos))
    pbcount <- 0

    for(i in gpos){
      pbcount <- pbcount + 1
      utils::setTxtProgressBar(pb, pbcount)
      bn <- decipherGaps_row(blast, i)
      blast <- rbind(blast, bn)
    }
    close(pb)
    blast <- blast[-gpos, ] #remove old rows

  } else{
    #parallel

    requireNamespace("foreach")
    requireNamespace("doParallel")

    `%dopar%` <- foreach::`%dopar%`

    ncores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, cores = ncores)
    print(paste("Paralellization registered with cores = ", ncores))


    bn <- foreach::foreach(i = gpos, .combine = rbind, .packages = c("magrittr")) %dopar% {
      decipherGaps_row(blast, i)
    }

    parallel::stopCluster(cl)

    blast <- rbind(blast, bn)
    blast <- blast[-gpos, ]

  }


}


decipherGaps_row <- function(blast, i){

  #run pairwise alignment
  bl.i <- blast[i,]
  qfrag <- with(bl.i,substr(qSeq,qStart,qEnd))
  sfrag <- with(bl.i,substr(sSeq,sStart,sEnd))

  msabl <- utils::capture.output(msa::msaClustalW(c(qfrag,sfrag) %>%
                                                    as.character %>% (Biostrings::AAStringSet)))
  msabl <- msabl[(1:2)+length(msabl)-3]
  msagap <- gsub("\\[\\d\\] ","",msabl)

  if(nchar(msagap[1]) != nchar(msagap[2])){
    stop("Error: decipherGaps: aligning two sequences with different nchar")}

  #make a note of gaps
  g1 <- gregexpr("-", msagap[1])
  g2 <- gregexpr("-", msagap[2])
  g <- c(unlist(g1), unlist(g2))

  #identify ungapped fragments
  s <- vector(); e <- vector() #keep track of fragment start/end

  if(1 %in% g){gap <- TRUE #check whether starts on a gap
  } else{s <- 1; gap <- FALSE}

  for(j in 2:nchar(msagap[1])){ #loop through remaining & make a note of switch
    if(j %in% g){
      #when a new gap starts, mark previous as end
      if(gap == FALSE){e <- c(e, j-1)}
      gap <- TRUE
    } else{
      #when a new non-gap starts, mark as start
      if(gap == TRUE){s <- c(s, j)}
      gap <- FALSE
    }
    if((j == nchar(msagap[1])) & (gap == FALSE)){e <- c(e, j)}
  }

  aln <- data.frame(cbind(s, e))
  aln.q <- aln.s <- data.table::setnames(data.frame(matrix(nrow=0, ncol=2)), names(aln))

  #map fragments to positions in blast alignment table
  for(j in 1:nrow(aln)){
    #need to subtract any gaps that precede an alignment fragment
    g.q <- gregexpr("-", substr(msagap[1], 1, aln[j, 1])) %>% unlist
    g.s <- gregexpr("-", substr(msagap[2], 1, aln[j, 1])) %>% unlist

    if(g.q[1] == -1){g.q <- 0} else{g.q %<>% length}
    if(g.s[1] == -1){g.s <- 0} else{g.s %<>% length}

    aln.q %<>% rbind(aln[j, ] - g.q)
    aln.s %<>% rbind(aln[j, ] - g.s)
  }

  #incrememnt based on alignment start position in bl.i
  aln.q <- aln.q + bl.i$qStart-1
  aln.s <- aln.s + bl.i$sStart-1

  #update blast alignment table
  bn <- blast[i, ]
  for(j in 1:nrow(aln)){
    bn[1, c("Gaps", "qStart", "qEnd", "sStart", "sEnd")] <-
      c(0, aln.q[j, ], aln.s[j, ]) %>% data.frame

    return(bn)
  }

}