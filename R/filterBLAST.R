#' filterBLAST
#' Remove singletons, sequneces with no non-self alignments,
#' from BLAST alignment table.
#' @param input BLAST alignment table from which to remove singletons.

filterBLAST <- function(input){ #remove singletons
  output.dir <- glParamGet("output.dir")
  epath <- paste0(output.dir, "epitopes/")
  spath <- paste0(epath, "singletons.csv")

  #remove self alignments
  filter <- input[(input$qID != input$sID) | (grepl("\\.", input$qID)), ]
  s <- setdiff(input, filter)
  s <- s[!(s$qID %in% filter$qID), ]
  if(!dir.exists(epath)) dir.create(epath)
  if(file.exists(spath)){
    sing <- fread(spath, data.table=FALSE)
    write.csv(rbind(s,sing) %>% unique, spath, row.names = FALSE)
  }
  write.csv(s, spath, row.names = FALSE)


  return(filter)
}
