#' organizeBLAST
#'
#' Rename unwieldy columns from rBLAST output table.
#'
#' @param input Input default BLAST table.

organizeBLAST <- function(input){
  names.orig <- c("QueryID", "SubjectID", "Gap.Openings",
                  "Q.start", "Q.end", "S.start", "S.end")
  names.rep <- c("qID", "sID", "Gaps",
                 "qStart", "qEnd", "sStart", "sEnd")

  #remove unecessary columns
  for(i in 1:length(names.orig)){names(input)[
    names(input) == names.orig[i]] <- names.rep[i]}


  names.rem <- c("Perc.Ident", "Alignment.Length", "Mismatches")
  for(i in 1:length(names.rem)){
    if(names.rem[i] %in% names(input)){
      input <- input[, -names.rem[i], with=FALSE]
    }
  }
  # if(ncol(input) == 12)
  # input <- input[, -c(3:5)]
  return(input)
}
