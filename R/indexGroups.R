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

        # prepare list of each subdivided group
        current.plus.aln <- list()
        for(j in 1:nrow(dup.combn)){
          current.plus.aln[[j]] <-dup.combn[j,] %>% as.matrix %>% as.character
        }

      } else{
        current.plus.aln <- list(current.plus.aln)
      }



      #if current.plus.aln is a list, iterate through each element. otherwise run once
      # if any of these are represented in at least one group, merge all
      # representing groups along with current.plus.aln

      for(j in 1:length(current.plus.aln)){
        # find indices of groups containing at least one of the members
        aligning.groups <- sapply(current.plus.aln[[j]], function(x){
          grep(x, groups, fixed = TRUE)}) %>% unlist %>% unique

        if(length(aligning.groups) > 0){
          merged.group <- c(current.plus.aln[[j]],
                            groups[aligning.groups] %>% unlist) %>%
            unique %>% sort

          groups[k] <- list(merged.group)
          groups <- groups[-aligning.groups]
        } else{
          groups[k] <- list(current.plus.aln[[j]])
        }
          k <- length(groups) + 1

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
