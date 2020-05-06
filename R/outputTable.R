#' Generates tables of reportable data on identified minimal alignments.
#'
#' @param blast BLAST alignment table to report.
#' @param fasta.initial Starting peptide sequences as AAStringSet.
#' @param groups Groups table of final epitopes to report.
#' @param msacs MSA consensus sequences for each group in order.
#' @param key_filename Epitope Key table filename to write.
#' @param summary_filename Epitope Summary table filename to write.
#'
#' @export

  outputTable <- function(blast, fasta.initial, groups,
                          msacs, key_filename, summary_filename){

    #get peptide names for which meaningful alignments were found by BLAST
    ep_peptides <- gsub("\\.[0-9]+\\.[0-9]+","",groups$ID) %>% unique

    #identify singleton peptides (peptides for which meaningful alignments were not found)
    sing <- fasta.initial[!(names(fasta.initial) %in% ep_peptides)]

    #merge singletons and aligning peptides into single table with ID, Seq, and Group Number
    if(length(sing) == 0){
      full <- groups
    } else{
      si <- data.frame(ID = paste(names(sing), 1, (as.character(sing) %>% nchar),
                                  sep = "."), Seq = as.character(sing),
                       Group = "NA")
      full <- rbind(groups, si)
    }
    rownames(full) <- c(1:nrow(full))


    #define and start populating output table
    cnames <- c("id", "start", "end", "ep_seq", "peptide_seq",
                "group_number", "group_name","group_size","msa_consensus")
    output <- data.frame(matrix("NA", nrow=nrow(full), ncol=length(cnames))) %>%
      data.table::setnames(cnames)


    output$id <- full$ID %>% gsub("\\.[0-9]+\\.[0-9]+$", "", .)
    output[, c("start", "end")] <-
      stringr::str_extract(full$ID, "[0-9]+\\.[0-9]+$") %>% strsplit("\\.") %>%
      unlist %>% matrix(nrow = 2) %>% t
    output[,"start"] %<>% as.numeric; output[,"end"] %<>% as.numeric
    output$ep_seq <- full$Seq %>% as.character

    output$peptide_seq <- fasta.initial[output$id] %>% as.character

    #Group Info: group number, group name, msa cosnensus sequence
    output$group_number <- full$Group
    gmax <- max(groups$Group)
    for(i in 1:gmax){
      gnames <- groups$ID[groups$Group == i]
      gname <- paste(gnames, collapse = "__")
      output$group_name[output$group_number == i] <- gname
      output$group_size[output$group_number == i] <- length(gnames)
    }

    output$msa_consensus <- suppressWarnings(msacs[output$group_number %>% as.numeric])


    # sort
    output <- output[(order(output$id,
                            as.numeric(output$start),
                            as.numeric(output$end))),]

    output$position <- paste0(output$start, "_", output$end)
    output$id_group <- paste0(output$id, "_", output$group_number)

    epitope_key <- unique(output[, c("group_number",
                                     "group_size",
                                     "msa_consensus")]) %>% (stats::na.omit)
    epitope_key <- epitope_key[order(as.numeric(epitope_key$group_number)),]

    epitope_summary <- output %>%
      stats::aggregate(position ~ id_group, data = ., FUN = function(x){
        paste(x, collapse = "_")})
    epitope_summary <- epitope_summary %>%
      dplyr::mutate(id = (epitope_summary$id_group %>% gsub("_[0-9]+$","", .) %>% gsub("_NA$","", .))) %>%
      dplyr::mutate(group_number = (
        stringr::str_extract(epitope_summary$id_group, "[0-9]+$")) %>% as.numeric) %>%
      # dplyr::select(output$id, output$position, output$group_number) %>%
      dplyr::select(id, position, group_number) %>%
      tidyr::spread(group_number, position)



    output_data <- list(epitope_key, epitope_summary)
    names(output_data) <- c("epitope_key", "epitope_summary")

    return(output_data)
  }
