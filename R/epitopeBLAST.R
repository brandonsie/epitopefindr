#' epitopeBLAST
#' Identifies next index peptide, calls indexEpitopes to parse alignments to
#' that index peptide, and updates BLAST table accordingly.
#' @param path Path to current fasta file of sequences to process.
#' @param verbose Logical of whether or not to print status updates to console.
#' @export

epitopeBLAST <- function(path, verbose = FALSE){
  # == == == == == A. Load fasta, separate full peptides & epitopes. == =
  peptides <- Biostrings::readAAStringSet(path)
  epitopes <- peptides[grepl("\\.", names(peptides))]

  blast.main <- glGet("blast.main")
  gl <- c("output.dir","proj.id")
  for(i in gl){assign(i,glParamGet(i))}


  #subset out prior epitopes. ep.genes = gene|tile names of epitopes.
  ep.names <- names(epitopes) %>% strsplit("__") %>% unlist %>% as.character
  ep.genes <- sapply(1:length(ep.names), function(x){
    ep.names[x] %>% strsplit("\\.") %>% unlist %>% `[[`(1)})
  blast.current <- blast.main[!(blast.main$qID %in% ep.genes), ]

  # == == == == == B. Choose index peptide & imsadentify its epitopes. == =
  index <- blast.current %>% chooseIndex()
  ipath <- paste0(output.dir,"epitopes/","indexOrder.txt")
  write(index,ipath,append=TRUE)

  blast.index <- rbind(
    blast.main[blast.main$qID==index,-"nAlign"],
    qsSwap(blast.main[blast.main$sID==index,-"nAlign"])) %>% unique
  blast.backup <- blast.main

  if(!exists("verbose")){verbose <- FALSE}
  if(verbose == FALSE){
    epIndex <- blast.index %>% indexEpitopes
  } else{
    print(path)
    print(paste(index, ", Blast Rows:", nrow(blast.index)))
    system.time(epIndex <- blast.index %>% indexEpitopes %>% print)
    print(as.character(epIndex[, 1]))
    cat("\n")
  }

  # == == == == == C. Write new epitopes .fasta file. == =
  #replace index peptides with index epitopes
  blast.remain <- blast.current[blast.current$qID!=index, ] #aln, q != index
  epFrag <- data.frame(ID = names(epitopes), Seq = epitopes %>% as.character)
  epList <- data.frame(ID = blast.remain$qID, Seq = blast.remain$qSeq) %>%
    rbind(epIndex, epFrag) %>% unique #blast.remain + old eps + new index eps
  epList %<>% mergeFastaDuplicates #remove duplicates, merge nomenclature

  #write alignments: index epitopes, remaining peptides, former epitopes
  fpath <- paste0(output.dir, "epitopes/")
  if(!dir.exists(fpath)){dir.create(fpath)}
  if(path == paste0(output.dir, "epitopes/", proj.id, ".fasta")){
    path <- paste0(output.dir, "epitopes/epitopes1.fasta")
  } else{
    gpath <- substr(path, 3, nchar(path)) %>% strsplit("/") %>% unlist
    j <- readr::parse_number(gpath[length(gpath)]) + 1
    path <- paste0(output.dir, "epitopes/epitopes", j, ".fasta")
  }
  writeFastaAA(epList, path)
  return(path)
} #END epitopeBLAST()
