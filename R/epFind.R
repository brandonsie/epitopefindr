#' epFind
#' All-in-one function to call the major steps of epitopefindr.
#' @param proj.id Name of peptide cohort. Output written to output/proj.id/.
#' @param e.thresh Maximum e value to allow from blast.alignments. Values
#' above 1 are not fully supported.
#' @param g.method Grouping method, either "any" or "all". See groupMSA
#' documentation for further details.
#' @param aln.size Shortest alignment to consider valid. Smaller alignments
#' are ignored and do not influence final epitopes.
#' @param autorun Logical whether to excute functions or just initialize BLAST.
#' @param relocate.input Logical whether to move input fasta file to new location.
#' @param relocate.dir Directory into which to move files if relocate.input == TRUE.
#' @export



epFind <- function(proj.id,
                   e.thresh = 0.01,
                   g.method = "any",
                   aln.size = 7,
                   autorun = TRUE,
                   relocate.input = FALSE,
                   relocate.dir) {
  # == == == == == Setup/configuration steps. == == == == ==
  options(stringsAsFactors = FALSE)
  if (grepl("\\.fasta", proj.id)) {
    proj.id <- gsub("\\.fasta", "", proj.id)
  }

  printSignpost(0, c(proj.id, e.thresh, g.method))

  printSignpost(1)
  epSetupDirectory(proj.id, e.thresh, g.method) #prepare output directories
  epSetupPeptides() #cleanup input sequences
  epcontinue <-
    epSetupBLAST() #blast input seqs against each other, tidy data

  if (!epcontinue) {
    dbCleanup(relocate.input = relocate.input, relocate.dir = relocate.dir)
    return(paste("No BLAST alignments identified for", proj.id))
  }

  # load some data from global environment
  blast.main <- glGet("blast.main")


  gl <- c("path", "blast.id3", "blast.id4", "g.method")
  for (i in gl) {
    assign(i, glParamGet(i))
  }

  # == == == == == Main script execution. == == == == ==
  if (autorun) {
    0

    printSignpost(2, blast.main$qID %>% unique %>% length)
    path %<>% pbCycleBLAST(ncycles = "max")
    data.table::fwrite(glGet("blast.main"), blast.id3)

    printSignpost(3)
    path %<>% trimEpitopes()
    data.table::fwrite(glGet("blast.main"), blast.id4)

    printSignpost(4)
    indexGroups(path, mode = g.method) #group alinging eps

    printSignpost(5)
    groupMSA() #generate multiple sequence alignment motifs

    printSignpost(6)
    outputTable() #generate output table
    outputFiles() #copy relevant output files to a new directory


    #tidy up R environment and temporary files
    dbCleanup(relocate.input = relocate.input, relocate.dir = relocate.dir)

    printSignpost(7)
    return(glParamGet("output.dir"))
  }
}
