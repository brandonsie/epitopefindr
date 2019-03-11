#' outputFiles
#'
#' Generates some output files to report from the epitope search.
#'
#' @export

outputFiles <- function(){
  gl <- c("output.dir","proj.id","e.thresh","g.method")
  for(i in gl){assign(i,glParamGet(i))}

  epath <- paste0(output.dir, "epitopes/", proj.id, ".fasta")
  file.copy(epath, paste0(output.dir,"initial_peptides.fasta"))

  p1 <- c("Project", "E Threshold", "Group Method")
  p2 <- c(proj.id, e.thresh, g.method)
  pc <- data.table::setnames(cbind(p1, p2) %>% data.frame, c("parameter", "value"))
  data.table::fwrite(pc, paste0(output.dir, "parameters.txt"))

}
