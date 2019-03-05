#' epPB
#'
#' Short for epitope progress bar, epPB calls utils::txtProgressBar with some
#' standardized parameters to reduce redundancy.
#'
#' @param min txtProgressBar minimum value
#' @param max txtProgressBar maximum value
#'
#' @importFrom utils txtProgressBar

epPB <- function(min,max){
  txtProgressBar(min = min, max = max, initial = min,
                 char = "=", width = NA, style = 3, file = "")
}
