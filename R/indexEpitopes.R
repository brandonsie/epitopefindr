#' indexEpitopes
#' For a specified index peptide, identify the maximal intervals that represent
#' the consensus overlap of reported BLAST alignments to that peptide. Update
#' the BLAST table accordingly.
#' @param blast.index Name of index peptide to process.
#' @export

indexEpitopes <- function(blast.index){
  # == == == == == A. Set up data frames. == =
  pos00 <- blast.index[, c("qStart", "qEnd")] %>% unique %>% data.frame #
  pos00 <- pos00[(pos00$qEnd - pos00$qStart > 5), ]
  pos00 <- pos00[order(pos00$qEnd, pos00$qStart), ]	#primary sort by qEnd!
  rownames(pos00) <- c(1:nrow(pos00))

  pos0 <- pos00
  pos1 <- data.table::setnames(data.frame(matrix(nrow = 0, ncol = 2)), names(pos00))
  df <- IRanges::findOverlaps(makeIR(pos00), minoverlap = 7) %>% data.frame
  prev <- matrix(nrow = nrow(pos00), ncol = 1)

  # == == == == == Run overlap identification cycle == =
  for (i in 1:nrow(pos0)) {
    #for each aln, starting with earliest end
    #get list of all overlaps, sorted by query then by start pos of subject

    if ((pos0$qEnd[i] - pos0$qStart[i] > 5)) {
      # == == == == == First Epitope == == == == ==
      if (i == 1) {
        #update sdf and psd. this must come after the .pre assignment
        sh0 <- sh.new <- df$subjectHits[df$queryHits == i] %>% sort #aln overlaps
        psd.new <- pos0[sh.new, ] #take corresponding start/end from pos0

        #produce minimal epitope
        minep <- c(max(psd.new$qStart), min(psd.new$qEnd))
        pos1 <- data.table::setnames(rbind(pos1, minep), names(pos00))

        # crop other alns that start before (earliest end-5) to start there
        minst <- min(psd.new$qEnd) - 5
        pos0$qStart[pos0$qStart < minst] <- minst
        pos0[pos0$qEnd - pos0$qStart <= 5, ] <- c(0, 0) #remove from cropped

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
          minst <- min(psd.new$qEnd) - 5
          pos0$qStart[pos0$qStart < minst] <- minst

          #remove alignments with cropped length < 7 (incl earliest end aln)
          pos0[pos0$qEnd - pos0$qStart <= 5, ] <- c(0, 0)

        }
      }
    }
  }

  # Tidy up new overlap table
  posN <- (IRanges::countOverlaps(makeIR(pos00), makeIR(pos1), minoverlap = 7)) == 0
  w1 <- paste("CAUTION:",
              "indexEpitopes/posN len > 0 for index", blast.index$qID[1])
  if(length(posN[posN == TRUE])>0) {print(w1)}
  pos1 <- rbind(pos1, pos00[posN, ]) %>% unique
  gpos <- pos1

  # == == == == == C. Update blast.main. == =
  #amend blast.main positions to be trimmed to corresponding overlap index ep

  blast.main <- glGet("blast.main")
  for(i in 1:nrow(blast.index)){ #for each
    if(gpos[gpos$qStart == blast.index$qStart[i] &
            gpos$qEnd == blast.index$qEnd[i],] %>% nrow == 0){
      #ignore exact matches

      qpos <- c(1:nrow(blast.main))[blast.main$qID == blast.index$qID[i] &
                                      blast.main$sID == blast.index$sID[i] &
                                      blast.main$qStart == blast.index$qStart[i] &
                                      blast.main$qEnd == blast.index$qEnd[i] &
                                      blast.main$sStart == blast.index$sStart[i] &
                                      blast.main$sEnd == blast.index$sEnd[i]]
      spos <- c(1:nrow(blast.main))[blast.main$qID == blast.index$sID[i] &
                                      blast.main$sID == blast.index$qID[i] &
                                      blast.main$qStart == blast.index$sStart[i] &
                                      blast.main$qEnd == blast.index$sEnd[i] &
                                      blast.main$sStart == blast.index$qStart[i] &
                                      blast.main$sEnd == blast.index$qEnd[i]]

      pos <- c(qpos, spos) %>% na.omit %>% as.numeric


      if(length(pos)>0){

        gOverlap <- isOverlapping(blast.index[i,c("qStart","qEnd")],gpos)
        for(j in gOverlap){ #for each position of gpos that overlaps

          dstart <- gpos[j, 1] - blast.index[i, "qStart"] #g-b. positive
          dend <- gpos[j, 2] - blast.index[i, "qEnd"] #g-b. negative


          #update so that multiple rows are made for each j that aligns
          b.add <- blast.main[pos,]
          for(p in 1:length(pos)){
            b.add[p,"qStart"] %<>% + dstart
            b.add[p,"qEnd"] %<>% + dend
            b.add[p,"sStart"] %<>% + dstart
            b.add[p,"sEnd"] %<>% + dend
          }
          blast.main <- rbind(blast.main, b.add)


        }

        #then remove original qpos spos rows
        blast.main <- blast.main[-pos, ] %>% unique

      }
    }
  }

  blast.main %<>% removeSmallAln
  blast.main <- blast.main[!duplicated(
    blast.main[,c("qID","sID","qStart","qEnd","sStart","sEnd")]),]
  blast.main %<>% numAlignments()
  glAssign("blast.main", blast.main)

  # == == == == == D. Convert to epitope sequences & output. == =
  indexep <- data.frame(ID=character(), Seq=character())
  for(i in 1:nrow(gpos)){
    ID <- paste(blast.index$qID[1] %>% as.character,
                gpos[i, 1], gpos[i, 2], sep=".")
    Seq <- blast.index$qSeq[1] %>% substr(gpos[i, 1], gpos[i, 2])
    indexep %<>% rbind(data.table::setnames(as.data.frame(t(c(ID, Seq))), names(indexep)))
  }

  return(indexep)

} #end indexEpitopes
