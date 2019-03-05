#' makeIR
#'
#' converts a two-column data frame to an IRanges object for overlap calculations.
#'
#' @param df Data frame to convert to IR object.


makeIR <- function(df){
  if(class(df)[1] == "numeric" & length(df) == 2){
    df <- data.frame(df) %>% t

  } else{df <- data.frame(df)}

  IRanges::IRanges(df[, 1], df[, 2])
}
