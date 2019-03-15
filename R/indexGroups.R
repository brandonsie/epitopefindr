#' indexGroups
#'
#' Group trimmed sequences into aligning groups
#'
#' @param blast BLAST alignment table to process.
#' @param fasta Epitope sequence AAStringset object to process.
#' @param mode Grouping method. "any" creates a smaller number of groups, with
#' individual groups tending to have more members, such that in each group,
#' each member must align with at least 1 "any" of the other members. "all"
#' creates a larger number of groups, with individual groups tending to have
#' fewer members, wuch taht in each group, each member must align with "all"
#' other members.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#'
#' @export

indexGroups <- function(blast, fasta, mode="any", aln.size){

  # == == == == == A. Load final epitope list and blast table(s). == =
  epitopes <- unmergeFastaDuplicates(fasta)
  blast <- blast[blast$qEnd-blast$qStart >= (aln.size - 1) &
                             blast$sEnd-blast$sStart >= (aln.size - 1), ]

  # == == == == == B. For each epitope, find aligning epitopes. == =

  #define data frame from blast that represents each epitope's alignments
  b.simp <- data.table::setnames(data.frame(matrix(ncol=2, nrow=nrow(blast))), c("q", "s"))
  b.simp$q <- sapply(1:nrow(blast), function(x){
    paste(blast[x, c("qID", "qStart", "qEnd")], collapse=".")})
  b.simp$s <- sapply(1:nrow(blast), function(x){
    paste(blast[x, c("sID", "sStart", "sEnd")], collapse=".")})
  b.simp <- b.simp[b.simp$q %in% names(epitopes) &
                     b.simp$s %in% names(epitopes), ] %>% unique

  # == == == == == C. Merge above lists either inclusive or exclusive == =

  if(mode == "all"){
    # ---------- "Align to All" = Exclusive Groups (more, tighter) ----------
    #all members of a group must align with all other group members

    #define matrix of which fragments align with which fragments via b.sip
    galign <- matrix(0, nrow=length(epitopes), ncol=length(epitopes))
    for(i in 1:nrow(galign)){galign[i, i] <- 1} #set identity to 1
    for(i in 1:nrow(b.simp)){ #populate each cell corresponding to a b.simp row
      r = (1:length(epitopes))[names(epitopes) == b.simp$q[i]]
      c = (1:length(epitopes))[names(epitopes) == b.simp$s[i]]
      galign[r, c] <- 1
    }
    rownames(galign) <- colnames(galign) <- names(epitopes)

    #for each peptide, define its minimal group(s)
    aln <- list(Alignments = character())
    addAln <- function(aln, toadd){
      if((aln[[1]] %>% length) == 0){aln[1] <- list(toadd)
      } else{aln[(length(aln)+1)] <- list(toadd)}
      return(aln)
    }

    for(i in 1:nrow(galign)){ #
      g1 <- galign[galign[i, ] == 1, galign[i, ] == 1]
      if(mean(g1)!=1){ #if there's misalignment within this group, subdivide
        for(j in 1:nrow(g1)){ #for each peptide in g1 group
          g2 <- j
          for(k in 1:nrow(g1)){ #check all elements, add if it keeps mean = 1
            if(g1[c(g2, k), c(g2, k)]%>%mean == 1){g2 %<>% c(k) %>% unique %>% sort}
          }
          # print(g2)
          aln %<>% addAln(list(rownames(g1)[g2]))
        }
      } else{aln %<>% addAln(list(names(epitopes)[galign[i, ] == 1]))}
    }

    #remove duplicate groups and subset groups. order by length
    aln %<>% unique
    len <- sapply(1:length(aln), function(x){aln[x] %>% unlist %>% length})
    aln <- aln[order(-len)]
    rem <- as.vector(matrix(0, nrow=1, ncol=length(aln)))

    for(i in 1:(length(aln)-1)){for(j in (i+1):length(aln)){ #long i, shorter j
      long <- aln[i] %>% unlist
      short <- aln[j] %>% unlist
      if((short %in% long) %>% mean == 1){
        # print(paste(i, j))
        rem[j] <- 1
      }
    }}
    aln <- aln[!rem]
    groups <- aln

  } #end "all"

  if(mode == "any"){
    # ---------- "Align to Any" = Inclusive Groups (fewer, looser) ----------
    # all members of a group must align with at least one other group member.
    #merge groups until above condition is satisfied

    #define data frame to keep track of each epitope's alignments
    #populate ep.align$Align with lists of aligning epitopes based on b.simp
    ep.align <- data.table::setnames(
      data.frame(matrix(ncol=2, nrow=length(epitopes))),c("Ep", "Align"))
    ep.align$Ep <- names(epitopes)
    for(i in 1:nrow(ep.align)){
      ep.align$Align[i] <-
        c(ep.align$Ep[i], b.simp$s[b.simp$q == ep.align$Ep[i]]) %>% list}

    k <- 0
    pb <- epPB(min = -(ep.align$Align %>% unique %>% unlist %>% length) - 1,
               max = -length(epitopes))
    while((ep.align$Align %>% unique %>% unlist %>%
           length > length(epitopes))){
      k %<>% +1; if(k > 100) stop("ERROR: stuck in a while loop.")
      setTxtProgressBar(pb, -(ep.align$Align %>% unique %>% unlist %>% length))

      for(i in 1:nrow(ep.align)){
        eg <- sapply(1:nrow(ep.align), function(x){
          if(ep.align$Ep[i] %in% unlist(ep.align$Align[x])){
            return(1)
          } else return(0)
        }) %>% grep(1, .)
        ep.align$Align[i] <- ep.align[eg, "Align"] %>% unlist %>%
          unique %>% sort %>% list
      }
    }
    setTxtProgressBar(pb, -length(epitopes))

    #take final unique groups, sorted with longest groups first
    groups <- ep.align$Align %>% unique
    len <- sapply(1:length(groups), function(x) groups[x] %>% unlist %>% length)
    groups <- groups[order(-len)]
  } #end "any"


  # == == == == == D. Write groups information to table == =

  output <- data.frame(ID = names(epitopes), Seq = as.character(epitopes),
                       Group = NA)
  for(i in 1:length(groups)){
    output$Group[output$ID %in% groups[[i]]] <- i
  }

  return(output)

  # #old writing individual group files to path
  # gpath <- "groups/"
  # if(!dir.exists(gpath)) dir.create(gpath)
  # for(i in 1:length(groups)){
  #   ep <- epitopes[names(epitopes) %in% (groups[i] %>% unlist)]
  #   ID <- names(ep); Seq <- as.character(ep)
  #   gdf <- data.frame(ID, Seq)
  #   writeFastaAA(gdf, paste0(gpath, "group", i, ".fasta"))
  # }
  # print(paste(length(groups), "groups identified."))
  # return(length(groups))
}
