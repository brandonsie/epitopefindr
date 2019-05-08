#' All-in-one function to call the major steps of epitopefindr.
#'
#' @param data Biostrings::AAStringset input sequences to search for epitopes, or path to corresponding .fasta file.
#' @param output.dir Directory to which output files should be written.
#' @param e.thresh Maximum e-value to consider from BLASTp alignments of 'data'.
#' @param g.method Grouping method of alignments. Either 'any' or 'all'. See ?indexGroups
#' @param aln.size Minimum length of alignment to consider from BLASTp alignments of 'data'.
#' @param min.groupsize Minimum number of peptides per group to require in order to print a group.
#' @param min.consensus.pos Minimum number of amino acid consensus positions required in order to print a group.
#' @param consensus.thresh Two decreasing numeric values of upper and lower thresholds for sequence consensus.
#' @param peptide.nchar Maximum of character from peptide name to use in msa output. Default 50. Starts from left.
#' @param msa.width Controls whether or not MSA images have fixed or dynamic width. By default, msa.width is set to "dynamic", which causes the document dimentions of the resultant image to be calculated based on the lentht of the peptide name and the number of amino acids in the sequence alignment. If msa.width is instead set to a numeric, then an MSA will be printed with a fixed with that number of inches. With 50-character peptide.nchar and a maximum expected sequence alignment of 45 positions, an msa.width of 12 is more than sufficient.
#' @param verbose Logical to print progress updates.
#' @param pdflatex Logical whether or not to produce PDF LaTeX figures using pdflatex
#' @param pdftk Logical whether or not to merge msa pdfs using staplr and pdftk
#' @param pdfuniter Logical whether or not to merge msa pdfs using staplr and pdfuniter
#' @param make.png Locial whether or not to convert PDF output to PNG.
#' @param name.msa Filename for output merged pdf of msa logos.
#' @param name.alignments Filename for output spreadsheet of peptide alignments.
#' @param name.epitopekey Filename for output spreadhseet of epitopes per peptide.
#' @param name.epitopesum Filename for output summary sheet of epitopes.
#' @param use.doParallel Logical whether or not to use doParallel parallelization.
#'
#' @export

epfind <- function(data = NULL, output.dir = NULL,
                    e.thresh = 0.01, g.method = "any", aln.size = 7,
                    min.groupsize = 2, min.consensus.pos = 1, consensus.thresh = c(75, 50),
                   peptide.nchar = 50, msa.width = "dynamic",
                    verbose = TRUE, pdflatex = TRUE, pdftk = TRUE, pdfuniter = FALSE, make.png = FALSE,
                    name.msa = "msa.pdf",
                    name.alignments = "finalAlignments.csv",
                    name.epitopekey = "epitopeKey.csv",
                    name.epitopesum = "epitopeSummary.csv",
                   use.doParallel = FALSE
                    ){

  # ----------------------------------------------------------------------------
  # Check parameters

  # throw error if data or output.dir are not defined
  if(is.null(data)){
    stop("Error: epfind param 'data' is undefined.")
  }

  if(is.null(output.dir)){
    stop("Error: epfind param 'output.dir' is undefined.")
  }

  # read 'data' if input path to .fasta file
  if(class(data)[1] == "character"){
    data <- Biostrings::readAAStringSet(data)
  }

  if(length(data) == 0){
    stop("Error: epfind: zero peptides input to data parameter.")
  }

  # setup directories
  temp.dir <- paste0(output.dir,"/intermediate_files/")
  if(!dir.exists(output.dir)){dir.create(output.dir)}
  if(!dir.exists(temp.dir)){dir.create(temp.dir)}

  options(stringsAsFactors = FALSE)
  # ----------------------------------------------------------------------------
  # Prepare sequences

  writeFastaAA(data, paste0(temp.dir,"fasta0.fasta"))


  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 1 of 6: Preparing BLAST alignment data from input sequences.",
        "\n")
  }

  fasta1 <- tidyFasta(data)
  f1.path <- paste0(temp.dir,"fasta1.fasta")
  writeFastaAA(fasta1, f1.path)

  peptide.name.map <- data.frame(original = names(data), modified = names(fasta1))
  data.table::fwrite(peptide.name.map,
                     paste0(temp.dir, "peptide_name_map.txt"), sep = "\t")

  blast1 <- selfBLASTaa(f1.path)
  b1.path <- paste0(temp.dir, "blast1.csv")
  data.table::fwrite(blast1, b1.path)
  if(nrow(blast1) == 0){
    stop("Error: epfind: no BLAST alignments found among input peptides.")
  }

  blast2 <- threshBLAST(blast1, e.thresh)
  b2.path <- paste0(temp.dir, "blast2.csv")
  data.table::fwrite(blast2, b2.path)
  if(nrow(blast2) == 0){
    stop(paste(
      "Error: epfind: no BLAST alignments with e value below:", e.thresh))
  }

  blast3 <- prepareBLAST(blast2, fasta1, aln.size, use.doParallel)
  b3.path <- paste0(temp.dir, "blast3.csv")
  data.table::fwrite(blast3, b3.path)
  if(nrow(blast3) == 0){
    stop(paste(
      "Error: epfind: no BLAST alignments after removing self-alignments."))
  }
  # ----------------------------------------------------------------------------
  # Process alignment overlaps

  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 2 of 6: Simplifying alignments to minimal number of overlapping intervals.",
        "\n")
  }

  blast4fasta <- pbCycleBLAST(blast3, fasta1, aln.size)
  blast4 <- blast4fasta[[1]]
  fasta4 <- blast4fasta[[2]]
  b4.path <- paste0(temp.dir, "blast4.csv")
  data.table::fwrite(blast4, b4.path)
  f4.path <- paste0(temp.dir, "fasta4.fasta")
  writeFastaAA(fasta4, f4.path)

  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 3 of 6: Trimming interval sequences.", "\n")
  }
  blast5fasta <- trimEpitopes(blast4fasta, aln.size)
  blast5 <- blast5fasta[[1]]
  fasta5 <- blast5fasta[[2]]
  b5.path <- paste0(temp.dir, "blast5.csv")
  data.table::fwrite(blast5, b5.path)
  f5.path <- paste0(temp.dir, "fasta5.fasta")
  writeFastaAA(fasta5, f5.path)

  blast6 <- renameBLAST(blast5, peptide.name.map)
  fasta6 <- renameFasta(fasta5, peptide.name.map)
  b6.path <- paste0(temp.dir, "blast6.csv")
  data.table::fwrite(blast6, b6.path)
  f6.path <- paste0(temp.dir, "fasta6.fasta")
  writeFastaAA(fasta6, f6.path)

  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 4 of 6: Grouping aligning sequences.", "\n")
  }
  groups <- indexGroups(blast6, fasta6, mode = g.method, aln.size)
  g.path <- paste0(temp.dir, "groups.csv")
  data.table::fwrite(groups, g.path)

  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 5 of 6: Generating multiple sequence alignment logos.", "\n")
  }

  m.path <- paste0(temp.dir, "msa/")
  if(!dir.exists(m.path)){dir.create(m.path)}
  g.path <- paste0(temp.dir, "groups.csv")
  groupMSA(groups, m.path, min.groupsize, min.consensus.pos,
           consensus.thresh, peptide.nchar, msa.width, pdflatex, pdftk, pdfuniter, make.png = make.png)

  if(verbose){
    cat("\n", "[", format(Sys.time(), "%R:%S"), "]",
        "Step 6 of 6: Preparing output files.", "\n")
  }

  if(pdftk | pdfuniter){file.copy(paste0(m.path,"msa.pdf"), paste0(output.dir,"/", name.msa))}
  file.copy(paste0(temp.dir, "blast6.csv"), paste0(output.dir, "/", name.alignments))

  msa.cs <- readLines(paste0(m.path,"consensusSequences.txt"))
  k.path <- paste0(output.dir,"/",name.epitopekey)
  s.path <- paste0(output.dir,"/",name.epitopesum)
  outputTable(blast6, data, groups, msa.cs, k.path, s.path)


}

#' @rdname epfind
#' @export
epFind2 <- epfind