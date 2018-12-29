#' cycleBLAST
#' Repeatedly call epitopeBLAST to iterate through all index peptides
#' until all peptides are converted to trimmed epitopes.
#' @param path Path to current fasta file.
#' @param pb txtProgressBar object to update
#' @param n Value to update txtProgressBar object.
#' @param verbose Logical whether or not to print status updates to console.
#' @export

cycleBLAST <- function(path, pb, n, verbose = FALSE){

  if(n > 0){ # n counts down to zero
    setTxtProgressBar(pb, -n) #update progress bar in console

    #Optionally print time at each iteration for debugging / speed optimization
    if(verbose == TRUE) {print(paste(Sys.time(), path, "start"))}

    #Check that there  are still unprocessed index peptides (name lacks period)
    peptides <- Biostrings::readAAStringSet(path)
    if(names(peptides) %>% grepl("\\.", .) %>% mean < 1){

      path <- epitopeBLAST(path, verbose)
      glParamAssign("path", path)
      path <- glParamGet("path")
      n <- n - 1 # -sum(!(grepl("\\.", names(epitopesI))))

      path <- cycleBLAST(path, pb, n)
    }
  }

  close(pb)

  return(path)
}
