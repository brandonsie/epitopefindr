#' dbCleanup
#'
#' Directory management for epitopefindr
#'
#' @param relocate.input Logical whether to move input fasta file to new location.
#' @param relocate.dir Directory into which to move files if relocate.input == TRUE.
#'
#' @export

dbCleanup <- function(relocate.input = FALSE, relocate.dir){

  output.dir <- glParamGet("output.dir")

  # # remove temporary files made during blast
  # epath <- paste0(output.dir, "epitopes/")
  # dbf <- list.files(epath)[list.files(epath) %>% grepl("psq|pin|phr", .)]
  # dbf <- paste0(output.dir, "epitopes/", dbf)
  # file.remove(dbf)

  pathi <- glParamGet("pathi")

  #if specified (off by default), move input fasta to a new directory
  if(relocate.input){
    if(!dir.exists(relocate.dir)){dir.create(relocate.dir)}
    file.copy(pathi,relocate.dir)
    file.remove(pathi)
  }

}
