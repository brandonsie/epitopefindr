#' printSignpost
#'
#' Some status updates that can be printed if verbose == TRUE.
#'
#' @param step.num Numerif reference of which step to display status for.
#' @param ... Extra parameters to past into step.num 0.
#'
#' @export

printSignpost <- function(step.num, ... = NA){
  if(step.num == 0){
    print(paste("Running epitope finder on [",...[1],".fasta] with [e < ",
                ...[2],"] and grouping method [",...[3],"].",sep=""))
  } else if(step.num == 1){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 1 of 6:",
                "BLASTing input sequences against each other."))
  } else if(step.num == 2){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 2 of 6:",
                "Identifying epitopes for", ...,"peptides."))
  } else if(step.num == 3){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 3 of 6:",
                "Looping back through peptides in reverse order."))
  } else if(step.num == 4){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 4 of 6:",
                "Grouping epitope sequences."))
  } else if(step.num == 5){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 5 of 6:",
                "Generating multiple sequence alignment motifs."))
  } else if(step.num == 6){
    print(paste(format(Sys.time(), "%H:%M:%S"),"Step 6 of 6:",
                "Preparing other output files."))
  } else if(step.num == 7){
    print(paste(format(Sys.time(), "%H:%M:%S"),"epitopeFinder run complete!",
                "Output files have been written to",glParamGet("output.dir")))
  } else {stop("Error: .printSignpost: invalid step number.")}
}

