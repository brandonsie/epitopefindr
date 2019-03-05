#' epSetupBLAST
#'
#' Run BLASTp on input peptides against each other & pre-process for use

epSetupBLAST <- function(){

  gl <- c("path","blast.id1","blast.id2","e.thresh")
  for(i in gl){assign(i,glParamGet(i))}

  #check for previously written blast table and only re-compute if not present
  if(file.exists(blast.id2)){
    print(paste("Loading existing blast table from", blast.id2))
    blast.main <- data.table::fread(blast.id2)

  } else {
    if(file.exists(blast.id1)){
      print(paste("Loading existing blast table from", blast.id1))
      blast.path <- data.table::fread(blast.id1)
    } else{
      print("Blasting starting sequences against each other...")
      blast.path <- data.table::data.table(selfBLASTaa(path)) #blast seq against eachother
      if(nrow(blast.path) == 0){return(FALSE)}
      data.table::fwrite(blast.path, blast.id1) #write blast to csv
    }

    blast.thresh <- blast.path[blast.path$E < as.numeric(e.thresh), ]
    if(nrow(blast.thresh) == 0){return(FALSE)}

    blast.main <- prepareBLAST(blast.thresh)
    if(nrow(blast.main) == 0){return(FALSE)}
    data.table::fwrite(blast.main, blast.id2)

  }

  glAssign("blast.main", blast.main)
  return(TRUE)
}
