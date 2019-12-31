#' Group trimmed sequences into aligning groups
#'
#' @param blast BLAST alignment table to process.
#' @param fasta Epitope sequence AAStringset object to process.
#' @param mode Grouping method. "any" creates a smaller number of groups, with
#' individual groups tending to have more members, such that in each group,
#' each member must align with at least 1 "any" of the other members. "all"
#' creates a larger number of groups, with individual groups tending to have
#' fewer members, such that in each group, each member must align with "all"
#' other members.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#' @param peptide.once.per.group Logical defining whether or not an epitope is constrained to only include one interval per contributing peptide. Default value TRUE. Only affects alignment mode = "any".
#'
#' @export

indexGroups <- function(blast, fasta, mode="any", aln.size, peptide.once.per.group = TRUE){

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


  if(mode == "any" & peptide.once.per.group == FALSE){
    # ---------- "Align to Any" = Inclusive Groups (fewer, looser) ----------
    # all members of a group must align with at least one other group member.
    #merge groups until above condition is satisfied
    # peptide.once.per.group == FALSE. allow one peptide to occur more than once in a given epitope

    groups <- list()
    k <- 1 #next group number

    pb <- epPB(min = 0, max = length(epitopes))

    for(i in 1:length(epitopes)){
      utils::setTxtProgressBar(pb, i)

      # Get current epitope and all epitopes that align to it
      current <- names(epitopes)[i]
      aln.current <- b.simp$s[b.simp$q == current]
      current.plus.aln <- c(current, aln.current) %>% unique %>% sort

      # if any of these are represented in at least one group, merge all
      # representing groups along with current.plus.aln

      # find indices of groups containing at least one of the members
      aligning.groups <- sapply(current.plus.aln, function(x){
        grep(x, groups, fixed = TRUE)}) %>% unlist %>% unique

      if(length(aligning.groups) > 0){
        merged.group <- c(current.plus.aln,
                          groups[aligning.groups] %>% unlist) %>%
          unique %>% sort



        groups[k] <- list(merged.group)
        groups <- groups[-aligning.groups]
      } else{
        groups[k] <-list(current.plus.aln)
      }

      k <- length(groups) + 1

    }
    close(pb)

    len <- sapply(1:length(groups), function(x) groups[x] %>% unlist %>% length)
    groups <- groups[order(-len)]


  } # end "any"


  if(mode == "any" & peptide.once.per.group == TRUE){
    # ---------- "Align to Any" = Inclusive Groups (fewer, looser) ----------
    # all members of a group must align with at least one other group member.
    #merge groups until above condition is satisfied
    # peptide.once.per.group == TRUE. allow one peptide to occur only once in a given epitope

    groups <- list()
    k <- 1 #next group number

    pb <- epPB(min = 0, max = length(epitopes))

    for(i in 1:length(epitopes)){
      utils::setTxtProgressBar(pb, i)
      print(paste("[[--new i--]]", i)) #(!) debug

      # Get current epitope and all epitopes that align to it
      current <- names(epitopes)[i]
      aln.current <- b.simp$s[b.simp$q == current]
      current.plus.aln <- c(current, aln.current) %>% unique %>% sort

      # (!) split by duplicate peptide
      pepnames <- current.plus.aln %>% strsplit("\\.") %>%
        lapply(FUN = function(x) x[[1]]) %>% unlist

      if(length(unique(pepnames)) < length(pepnames)){
        #if there's at least one duplicated peptide, subdivide groups
        duplicated.peps <- pepnames[duplicated(pepnames)] %>% unique
        dup.list <- list()
        for(j in 1:length(duplicated.peps)){
          dup.list[[j]] <- current.plus.aln[pepnames == duplicated.peps[j]]
        }
        dup.combn <- dup.list %>% expand.grid

        non.duplicated.peps <- current.plus.aln[!(pepnames %in% duplicated.peps)]

        # prepare list of each subdivided group
        current.plus.aln <- list()
        for(j in 1:nrow(dup.combn)){
          current.plus.aln[[j]] <-dup.combn[j,] %>% as.matrix %>% as.character %>% c(non.duplicated.peps)
        }

      } else{
        current.plus.aln <- list(current.plus.aln)
      }



      #if current.plus.aln is a list, iterate through each element. otherwise run once
      # if any of these are represented in at least one group, merge all
      # representing groups along with current.plus.aln
      #(!) need to prevent duplicate pep entry here

      print(paste("[length j]", length(current.plus.aln)))#(!) debug
      for(j in 1:length(current.plus.aln)){
        if(j%%20==0) print(paste("j", j)) #(!) debug
        # get peptide names of previoulsy grouped peptides for upcoming comparison
        if(length(groups) > 0){
          groups.pepnames <- lapply(groups, FUN = function(x){
            x %>% strsplit("\\.") %>% lapply(FUN = function(x) x[[1]]) %>% unlist
          }) #still returns a list unless unlist here

        } else{groups.pepnames <- ""}

        # find indices of groups containing at least one of the members and no duplicate peptides
        aligning.groups <- sapply(current.plus.aln[[j]], function(x){
          all.aln.groups <- sapply(x, function(x){
            grep(x, groups, fixed = TRUE)
          }) %>% unlist %>% unique
          #grep(x, groups, fixed = TRUE)

          this.pepname <- x %>% strsplit("\\.") %>% lapply(
            FUN = function(x) x[[1]]) %>% unlist

          pep.dup.groups <- sapply(
            this.pepname, function(x){grep(x, groups.pepnames, fixed = TRUE)}) %>%
            unlist %>% unique
          #grep(this.pepname, groups.pepnames, fixed = TRUE)

          return(all.aln.groups[!(all.aln.groups %in% pep.dup.groups)])
        }) %>% unlist %>% unique

        # if there are valid aligning groups, combine existing group with peptides aligning to current
        if(length(aligning.groups) > 0){
          for(m in 1:length(aligning.groups)){
            merged.group <- c(current.plus.aln[[j]],
                              groups[[aligning.groups[m]]]) %>% unique %>% sort
            groups[[aligning.groups[m]]] <- merged.group
          }
        } else{ #if there are no aligning groups, create new group
          groups[k] <- list(current.plus.aln[[j]])
          k <- k + 1
        }

      }



    }
    close(pb)

    len <- sapply(1:length(groups), function(x) groups[x] %>% unlist %>% length)
    groups <- groups[order(-len)]


  } # end any, true

  # == == == == == D. Write groups information to table == =

  output <- data.frame(ID = names(epitopes), Seq = as.character(epitopes),
                       Group = NA)
  for(i in length(groups):1){
    output$Group[output$ID %in% groups[[i]]] <- i
  }

  if(sum(is.na(output) > 0)){
    warning("Warning (epitopefindr::indexGroups) some aligning peptides not assigned to groups.")
    output <- na.omit(output)
  }

  return(output)

}
