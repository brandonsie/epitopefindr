#' unmergeFastaDuplicates
#'
#' If a fasta entry has a concatenated name, e.g. from a previous call to
#' mergeFastaDuplicates, then separate the concatenanted name and make
#' a separate entry for each.
#'
#' @param input XStringSet object with concatenated names.
#' @param sep Delimiter to search for and parse.
#' @return XStringSet object modified so that a single entries with delimited
#' names are replaced with multiple entries with identical sequences and
#' unmerged names.
#'
#' @export

unmergeFastaDuplicates <- function(input, sep = "__"){
  np <- input %>% names %>% stringr::str_count(sep) #number of peptides - 1
  for(k in 1:length(np)){input <- c(input, rep(input[k], np[k]))}
  for(k in 1:length(input)){	#rename duplicated peptides
    names(input)[names(input) == names(input)[k]] <-
      unlist(strsplit(names(input)[k], sep))}
  return(input)
}
