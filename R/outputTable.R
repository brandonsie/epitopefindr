#' outputTable
#' Generates spreadsheet of reportable data on identified minimal alignments.
#' @export

outputTable <- function(){
  gl <- c("output.dir","proj.id","path","path0","blast.id4","e.thresh","g.method")
  for(i in gl){assign(i,glParamGet(i))}


  #load epitope info
  epitopes <- Biostrings::readAAStringSet(path) %>% unmergeFastaDuplicates
  ep <- data.frame(ID = names(epitopes), Seq = as.character(epitopes))
  original <- Biostrings::readAAStringSet(path0)

  ## load singletons
  spath <- paste0(output.dir, "epitopes/singletons.csv")
  sing <- data.table::fread(spath)
  sf <- data.frame(ID = paste(sing$qID, sing$qStart, sing$qEnd, sep="."),
                   Seq = sing$qSeq)

  #merge epitopes and singletons
  full <- rbind(ep, sf); rownames(full) <- c(1:nrow(full))
  # haspipe <- ifelse(length(grep("\\|",full$ID)) > 0,TRUE,FALSE)

  #load msa consensus sequences
  cpath <- paste0(output.dir,"msa/consensusSequences.txt")
  cseqs <- gsub("-","",readLines(cpath))

  #define and start populating output table
  cnames <- c("gene", "tile", "start", "end", "seq_ep", "seq_tile",
              "grpnum", "grpname","grpsize","msaconsensus")
  output <- data.table::setnames(data.frame(matrix("NA", nrow=nrow(full),
                                       ncol=length(cnames))), cnames)

  for(i in 1:nrow(output)){
    if(full$ID[i] %>% grepl("\\|",.)){
      output[i,c("gene","tile")] <- strsplit(full$ID[i],"\\|") %>% unlist
    } else{output[i,c("gene","tile")] <- full$ID[i]}
  }

  output[, c("tile", "start", "end")] <- strsplit(output$tile, "\\.") %>%
    unlist %>% as.vector %>% matrix(nrow=3) %>% t
  output[,"start"] %<>% as.numeric; output[,"end"] %<>% as.numeric
  output$seq_ep <- full$Seq %>% as.character

  for(i in 1:nrow(output)){
    output$seq_tile[i] <- original[
      strsplit(full$ID[i],"\\.") %>% unlist %>% magrittr::extract(1)] %>%
      as.character
  }

  #load group info
  gpath <- paste0(output.dir, "groups/")
  gmax <- list.files(gpath) %>% (readr::parse_number) %>% max
  for(i in 1:gmax){
    gi <- Biostrings::readAAStringSet(paste0(gpath, "group", i, ".fasta"))
    gn <- paste(names(gi), collapse="__")

    output$grpnum[full$ID %in% names(gi)] %<>% paste(i, sep=", ")
    output$grpsize[full$ID %in% names(gi)] %<>% paste(length(gi), sep=", ")
    output$grpname[full$ID %in% names(gi)] <- gn
  }
  output$grpnum <- gsub("NA, ", "", output$grpnum)
  output$grpsize <- gsub("NA, ", "", output$grpsize)

  output$msaconsensus <- suppressWarnings(cseqs[output$grpnum %>% as.numeric])


  # sort
  output <- output[(order(output$gene,output$tile,
                          as.numeric(output$start),
                          as.numeric(output$end))),]

  opath <- paste0(output.dir, proj.id,"_e-",e.thresh,"_g-",g.method,
                  "_outputTable.csv")
  fwrite(output, opath)

}
