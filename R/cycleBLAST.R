#' cycleBLAST
#' Repeatedly call epitopeBLAST to iterate through all index peptides
#' until all peptides are converted to trimmed epitopes.
#' @param data List of BLAST alignment table and fastra file to process.
#' @param pb txtProgressBar object to update
#' @param n Value to update txtProgressBar object.
#' @param verbose Logical whether or not to print status updates to console.
#' @return List data containing modified BLAST alignment table and fasta file.
#' @export



cycleBLAST <- function(data, pb, n, verbose = FALSE){

  options(stringsAsFactors = FALSE)
  if(n > 0){ # n counts down to zero
    utils::setTxtProgressBar(pb, -n) #update progress bar in console

    #Optionally print time at each iteration for debugging / speed optimization
    if(verbose == TRUE) {print(paste(Sys.time(), path, "start"))}

    # #Check that there  are still unprocessed index peptides (name lacks period)
    # print(n)
    # print(length(names(data[[2]])))
    # print(names(data[[2]]) %>% grepl("\\.",.)) %>% mean

    if(names(data[[2]]) %>% grepl("\\.", .) %>% mean < 1){

      epData <- epitopeBLAST(data)
      n <- n - 1

      cycData <- cycleBLAST(epData, pb, n)
    }
  } else{
    return(data)
  }

}
