#' For a specified index peptide, identify the maximal intervals that represent
#' the consensus overlap of reported BLAST alignments to that peptide. Update
#' the BLAST table accordingly.
#'
#' @param blast BLAST alignment table to process.
#' @param index Name of index peptide to process.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#'
#' @export

indexEpitopes <- function(blast, index, aln.size){


  blast.index <- rbind(
    blast[blast$qID==index,-"nAlign"],
    qsSwap(blast[blast$sID==index,-"nAlign"])) %>% unique

  # == == == == == A. Set up data frames. == =
  pos00 <- blast.index[, c("qStart", "qEnd")] %>% unique %>% data.frame #
  pos00 <- pos00[(pos00$qEnd - pos00$qStart >= (aln.size - 1)), ]
  pos00 <- pos00[order(pos00$qEnd, pos00$qStart), ]	#primary sort by qEnd!
  rownames(pos00) <- c(1:nrow(pos00))

  pos0 <- pos00
  pos1 <- data.table::setnames(data.frame(matrix(nrow = 0, ncol = 2)), names(pos00))
  df <- IRanges::findOverlaps(makeIR(pos00), minoverlap = aln.size) %>% data.frame
  prev <- matrix(nrow = nrow(pos00), ncol = 1)

  # == == == == == Run overlap identification cycle == =
  for (i in 1:nrow(pos0)) {
    #for each aln, starting with earliest end
    #get list of all overlaps, sorted by query then by start pos of subject

    if ((pos0$qEnd[i] - pos0$qStart[i] >= (aln.size - 1))) {
      # == == == == == First Epitope == == == == ==
      if (i == 1) {
        #update sdf and psd. this must come after the .pre assignment
        sh0 <- sh.new <- df$subjectHits[df$queryHits == i] %>% sort #aln overlaps
        psd.new <- pos0[sh.new, ] #take corresponding start/end from pos0

        #produce minimal epitope
        minep <- c(max(psd.new$qStart), min(psd.new$qEnd))
        pos1 <- data.table::setnames(rbind(pos1, minep), names(pos00))

        # crop other alns that start before (earliest end-5) to start there
        minst <- min(psd.new$qEnd) - (aln.size - 2)
        pos0$qStart[pos0$qStart < minst] <- minst
        pos0[pos0$qEnd - pos0$qStart < (aln.size - 2), ] <- c(0, 0) #remove from cropped

      } else {
        # == == == == == Subsequent Epitopes == == == == ==
        prev[i - 1] <- sh0 %>% list #keep previous data
        sh0 <- sh <- df$subjectHits[df$queryHits == i] %>% sort #new aln overlp
        for(j in sh0[sh0<i]){sh <- sh[!(sh %in% unlist(prev[j]))]} #excl prev

        if(length(sh0[sh0<i]) == 0){sh <- c(sh, i) %>% unique %>% sort}
        sh.new <- sh

        #take alignments to i that didn't align to previous
        # sh.new <- rem[!(rem %in% c(1:nrow(pos0))[pos0$qEnd == 0])]
        if(length(sh.new)>0){
          sh.new <- c(sh.new, i) %>% unique %>% sort
          psd.new <- pos00[sh.new, ]

          # == == == == == produce minimal epitope == == == == ==
          minep <- c(max(psd.new$qStart), min(psd.new$qEnd))
          pos1 <- data.table::setnames(rbind(pos1, minep), names(pos00))

          # == == == == == crop alns that start before (end-5) == == == == ==
          minst <- min(psd.new$qEnd) - (aln.size - 2)
          pos0$qStart[pos0$qStart < minst] <- minst

          #remove alignments with cropped length < 7 (incl earliest end aln)
          pos0[pos0$qEnd - pos0$qStart <= (aln.size - 2), ] <- c(0, 0)

        }
      }
    }
  }

  # Tidy up new overlap table
  posN <- (IRanges::countOverlaps(makeIR(pos00), makeIR(pos1), minoverlap = aln.size)) == 0
  w1 <- paste("CAUTION:",
              "indexEpitopes/posN len > 0 for index", blast.index$qID[1])
  if(length(posN[posN == TRUE])>0) {print(w1)}
  pos1 <- rbind(pos1, pos00[posN, ]) %>% unique
  gpos <- pos1

  # == == == == == C. Update blast. == =
  #amend blast positions to be trimmed to corresponding overlap index ep

  for(i in 1:nrow(blast.index)){ #for each
    if(gpos[gpos$qStart == blast.index$qStart[i] &
            gpos$qEnd == blast.index$qEnd[i],] %>% nrow == 0){
      #ignore exact matches

      qpos <- c(1:nrow(blast))[blast$qID == blast.index$qID[i] &
                                      blast$sID == blast.index$sID[i] &
                                      blast$qStart == blast.index$qStart[i] &
                                      blast$qEnd == blast.index$qEnd[i] &
                                      blast$sStart == blast.index$sStart[i] &
                                      blast$sEnd == blast.index$sEnd[i]]
      spos <- c(1:nrow(blast))[blast$qID == blast.index$sID[i] &
                                      blast$sID == blast.index$qID[i] &
                                      blast$qStart == blast.index$sStart[i] &
                                      blast$qEnd == blast.index$sEnd[i] &
                                      blast$sStart == blast.index$qStart[i] &
                                      blast$sEnd == blast.index$qEnd[i]]

      pos <- c(qpos, spos) %>% (stats::na.omit) %>% as.numeric


      if(length(pos)>0){

        gOverlap <- isOverlapping(blast.index[i,c("qStart","qEnd")],gpos, aln.size)
        for(j in gOverlap){ #for each position of gpos that overlaps

          dstart <- gpos[j, 1] - blast.index[i, "qStart"] #g-b. positive
          dend <- gpos[j, 2] - blast.index[i, "qEnd"] #g-b. negative


          #update so that multiple rows are made for each j that aligns
          b.add <- blast[pos,]
          for(p in 1:length(pos)){
            b.add[p,"qStart"] %<>% + dstart
            b.add[p,"qEnd"] %<>% + dend
            b.add[p,"sStart"] %<>% + dstart
            b.add[p,"sEnd"] %<>% + dend
          }
          blast <- rbind(blast, b.add)


        }

        #then remove original qpos spos rows
        blast <- blast[-pos, ] %>% unique

      }
    }
  }

  blast %<>% removeSmallAln(aln.size)
  blast <- blast[!duplicated(
    blast[,c("qID","sID","qStart","qEnd","sStart","sEnd")]),]
  blast %<>% numAlignments()

  # == == == == == D. Convert to epitope sequences & output. == =
  indexep <- data.frame(ID=character(), Seq=character())
  for(i in 1:nrow(gpos)){
    ID <- paste(blast.index$qID[1] %>% as.character,
                gpos[i, 1], gpos[i, 2], sep=".")
    Seq <- blast.index$qSeq[1] %>% substr(gpos[i, 1], gpos[i, 2])
    indexep %<>% rbind(data.table::setnames(as.data.frame(t(c(ID, Seq))), names(indexep))) %>% as.data.frame
  }

  data <- list(blast = blast, indexep = indexep)
  return(data)

} #end indexEpitopes
