#' filterBLAST
#' Remove singletons, sequneces with no non-self alignments,
#' from BLAST alignment table.
#' @param blast BLAST alignment table from which to remove singletons.
#' @export

filterBLAST <- function(blast){ #remove singletons
  # if(is.null(output.dir)){
  #   spath <- "singletons.csv"
  # } else{spath <- paste0(output.dir,"singletons.csv")}

  #remove self alignments
  filter <- blast[(blast$qID != blast$sID) | (grepl("\\.", blast$qID)), ]
  # s <- setdiff(blast, filter)
  # s <- s[!(s$qID %in% filter$qID), ]
  # if(file.exists(spath)){
  #   sing <- fread(spath, data.table=FALSE)
  #   utils::write.csv(rbind(s,sing) %>% unique, spath, row.names = FALSE)
  # }
  # utils::write.csv(s, spath, row.names = FALSE)


  return(filter)
}
