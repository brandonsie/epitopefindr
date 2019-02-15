#' outputTable
#' Generates spreadsheet of reportable data on identified minimal alignments.
#' @param blast BLAST alignment table to report.
#' @param fasta.initial Starting peptide sequences as AAStringSet.
#' @param groups Groups table of final epitopes to report.
#' @param msacs MSA consensus sequences for each group in order.
#' @param filename Output table filename to write.
#' @export

outputTable <- function(blast, fasta.initial, groups,
                        msacs, filename){


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
  cnames <- c("gene", "tile", "start", "end", "ep_seq", "tile_seq",
              "group_number", "group_name","group_size","msaconsensus")
  output <- data.frame(matrix("NA", nrow=nrow(full), ncol=length(cnames))) %>%
    data.table::setnames(cnames)

  for(i in 1:nrow(output)){
    if(full$ID[i] %>% grepl("\\|",.)){
      output[i,c("gene","tile")] <- strsplit(full$ID[i],"\\|") %>% unlist
    } else{output[i,c("gene","tile")] <- full$ID[i]}
  }

  output[, c("tile", "start", "end")] <- strsplit(output$tile, "\\.") %>%
    unlist %>% as.vector %>% matrix(nrow=3) %>% t
  output[,"start"] %<>% as.numeric; output[,"end"] %<>% as.numeric
  output$ep_seq <- full$Seq %>% as.character

  for(i in 1:nrow(output)){
    output$tile_seq[i] <- original[
      strsplit(full$ID[i],"\\.") %>% unlist %>% magrittr::extract(1)] %>%
      as.character
  }

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

  output$msaconsensus <- suppressWarnings(msacs[output$group_number %>% as.numeric])


  # sort
  output <- output[(order(output$gene,output$tile,
                          as.numeric(output$start),
                          as.numeric(output$end))),]

  data.table::fwrite(output, filename)

}
