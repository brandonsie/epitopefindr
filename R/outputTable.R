#' Generates spreadsheet of reportable data on identified minimal alignments.
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


  #load epitope info
  original <- fasta.initial

  #get peptide names for which meaningful alignments were found
  ep_peptides <- gsub("\\.[0-9]+\\.[0-9]+","",groups$ID) %>% unique

  #(!) bookmark setup data frames
  sing <- original[!(names(original) %in% ep_peptides)]
  if(length(sing) == 0){
    full <- groups
  } else{
    si <- data.frame(ID = paste(names(sing), 1, (as.character(sing) %>% nchar),
                                sep = "."), Seq = as.character(sing),
                     Group = "NA")
    full <- rbind(groups, si)
  }
  rownames(full) <- c(1:nrow(full))


  #merge epitopes and singletons
  # haspipe <- ifelse(length(grep("\\|",full$ID)) > 0,TRUE,FALSE)


  #define and start populating output table
  cnames <- c("id", "start", "end", "ep_seq", "peptide_seq",
              "group_number", "group_name","group_size","msa_consensus")
  output <- data.frame(matrix("NA", nrow=nrow(full), ncol=length(cnames))) %>%
    data.table::setnames(cnames)

  # old when gene and tile were separate columns. now just id
  # for(i in 1:nrow(output)){
  #   if(full$ID[i] %>% grepl("\\|",.)){
  #     output[i,c("gene","tile")] <- strsplit(full$ID[i],"\\|") %>% unlist
  #   } else{output[i,c("gene","tile")] <- full$ID[i]}
  # }
#
#   output[, c("tile", "start", "end")] <- strsplit(output$tile, "\\.") %>%
#     unlist %>% as.vector %>% matrix(nrow=3) %>% t

  output$id <- full$ID %>% gsub("\\.[0-9]+\\.[0-9]+$", "", .)
  output[, c("start", "end")] <-
    stringr::str_extract(full$ID, "[0-9]+\\.[0-9]+$") %>% strsplit("\\.") %>%
    unlist %>% matrix(nrow = 2) %>% t
  output[,"start"] %<>% as.numeric; output[,"end"] %<>% as.numeric
  output$ep_seq <- full$Seq %>% as.character

  output$peptide_seq <- original[output$id] %>% as.character
  # for(i in 1:nrow(output)){
    # output$peptide_seq[i] <- original[
    #   strsplit(full$ID[i],"\\.") %>% unlist %>% magrittr::extract(1)] %>%
    #   as.character
  # }

  #Group Info: group number, group name, msa cosnensus sequence
  output$group_number <- full$Group
  gmax <- max(groups$Group)
  for(i in 1:gmax){
    gnames <- groups$ID[groups$Group == i]
    gname <- paste(gnames, collapse = "__")
    output$group_name[output$group_number == i] <- gname
    output$group_size[output$group_number == i] <- length(gnames)
  }

  # #(!) output table lost support for "all" grouping here
  # #load group info
  # gpath <- paste0(output.dir, "groups/")
  # gmax <- list.files(gpath) %>% (readr::parse_number) %>% max
  # for(i in 1:gmax){
  #   gi <- Biostrings::readAAStringSet(paste0(gpath, "group", i, ".fasta"))
  #   gn <- paste(names(gi), collapse="__")
  #
  #   output$group_number[full$ID %in% names(gi)] %<>% paste(i, sep=", ")
  #   output$group_size[full$ID %in% names(gi)] %<>% paste(length(gi), sep=", ")
  #   output$group_name[full$ID %in% names(gi)] <- gn
  # }
  # output$group_number <- gsub("NA, ", "", output$group_number)
  # output$group_size <- gsub("NA, ", "", output$group_size)

  output$msa_consensus <- suppressWarnings(msacs[output$group_number %>% as.numeric])


  # sort
  output <- output[(order(output$id,
                          as.numeric(output$start),
                          as.numeric(output$end))),]
  # output <- output[(order(output$gene,output$tile,
  #                         as.numeric(output$start),
  #                         as.numeric(output$end))),]

  #old output format
  # data.table::fwrite(output, "outputTable.csv")

  epitope_key <- unique(output[, c("group_number",
                                   "group_size",
                                   "msa_consensus")]) %>% na.omit
  epitope_key <- epitope_key[order(as.numeric(epitope_key$group_number)),]

  epitope_summary <- output %>%
    dplyr::mutate(position = paste0(start, "_", end)) %>%
    dplyr::mutate(id_group = paste0(id, "_", group_number)) %>%
    stats::aggregate(position ~ id_group, data = ., FUN = function(x){
      paste(x, collapse = "_")}) %>%
    dplyr::mutate(id = (id_group %>% gsub("_[0-9]+$","", .) %>% gsub("_NA$","", .))) %>%
    dplyr::mutate(group_number = (
      stringr::str_extract(id_group, "[0-9]+$"))) %>%
    dplyr::select(id, position, group_number) %>%
    tidyr::spread(group_number, position)

  # # re-sort summary
  # numeric.positions <- suppressWarnings(!(as.numeric(names(epitope_summary)) %>% is.na))
  # group.numbers <- names(epitope_summary)[numeric.positions] %>% as.numeric
  # to_re_sort <- epitope_summary[,numeric.positions]
  # to_re_sort <- to_re_sort[,order(group.numbers)]
  # epitope_summary[,numeric.positions] <- to_re_sort

  data.table::fwrite(epitope_key, key_filename)
  data.table::fwrite(epitope_summary, summary_filename)

}

# output$id <- full$ID %>% gsub("\\.[0-9]+\\.[0-9]+$", "", .)
# stringr::str_extract(full$ID, "[0-9]+\\.[0-9]+$") %>% strsplit("\\.") %>%
#   unlist %>% matrix(nrow = 2) %>% t
