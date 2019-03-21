#' Repeatedly call epitopeBLAST to iterate through all index peptides
#' until all peptides are converted to trimmed epitopes.
#'
#' @param data List of BLAST alignment table and fastra file to process.
#' @param pb txtProgressBar object to update
#' @param n Value to update txtProgressBar object.
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#' @param verbose Logical whether or not to print status updates to console.
#'
#' @return List data containing modified BLAST alignment table and fasta file.
#'
#' @export



cycleBLAST <- function(data, pb, n, aln.size, verbose = FALSE){

  options(stringsAsFactors = FALSE)

  utils::setTxtProgressBar(pb, -n) #update progress bar in console
  if(n > 0){ # n counts down to zero
    #Optionally print time at each iteration for debugging / speed optimization
    if(verbose == TRUE) {print(paste(Sys.time(), path, "start"))}

    # #Check that there  are still unprocessed index peptides (name lacks period)
    # print(n)
    # print(length(names(data[[2]])))
    # print(names(data[[2]]) %>% grepl("\\.",.)) %>% mean

    if(names(data[[2]]) %>% grepl("\\.", .) %>% mean < 1){

      data <- epitopeBLAST(data, aln.size)
      n <- n - 1

      data <- cycleBLAST(data, pb, n, aln.size)
    }
  }
  return(data)
}
