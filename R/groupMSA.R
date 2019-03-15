#' groupMSA
#'
#' Generate pdf multiple sequence alignment logos for identified aligning groups.
#'
#' @param groups Table of final peptides and pre-calculated groupings.
#' @param mpath Directory to write sequence alignment files.
#' @param min.groupsize Minimum number of peptides per group to require in order to print a group.
#' @param min.consensus.pos Minimum number of amino acid consensus positions required in order to print a group.
#' @param consensus.thresh Two decreasing numeric values of upper and lower thresholds for sequence consensus.
#' @param pdflatex Logical whether or not to produce PDF LaTeX figures using pdflatex
#' @param pdftk Logical whether or not to use staplr and pdftk to merge individual msa pdfs.
#' @param trim.groups Logical whether or not to apply msaTrim to edges of logos. Not implemented.
#' @param make.png Depreciated. Locial whether or not to convert PDF output to PNG.
#'
#' @export

groupMSA <- function(groups, mpath = "intermediate_files/msa/",
                     min.groupsize = 2, min.consensus.pos = 1, consensus.thresh = c(75, 50),
                     pdflatex = TRUE, pdftk = TRUE,
                     trim.groups = FALSE, make.png = FALSE){


  #setup output paths and directories
  cpath <- paste0(mpath,"consensusSequences.txt")
  if(!dir.exists(mpath)) {dir.create(mpath)}
  if(file.exists(cpath)) {file.remove(cpath)}

  #loop through each group
  num <- max(groups$Group)
  if(num < 1){stop("Error: no input groups specified.")}
  num.length <- nchar(num)

  pb <- epPB(0, num)

  for(i in 1:num){
    utils::setTxtProgressBar(pb, i)

    group <- Biostrings::AAStringSet(groups$Seq[groups$Group == i])
    names(group) <- groups$ID[groups$Group == i]
    group %<>% unmergeFastaDuplicates

    # == == == == == normalize peptide name length == == == == ==
    k = 50 #characters from peptide name to use

    #separate basenames and start/end positions
    basenames <- names(group) %>% gsub("\\.[0-9]+\\.[0-9]+$", "", .)
    positions <- stringr::str_extract(names(group), "\\.[0-9]+\\.[0-9]+$")

    #shorten long names, pad short names
    shortened.basenames <- substr(basenames, 1, k)
    lengthened.basenames <- stringr::str_pad(shortened.basenames,
                                             width=k+1,side="right",pad=".")

    #reappend start & end positions, update names in `group` object
    fully.modified.names <- paste0(lengthened.basenames, positions)
    names(group) <- fully.modified.names

    #trim group if too many members to msaprettyprint(). Store consensus sequence.
    group.vector <- as.character(group)
    mg <- Biostrings::AAMultipleAlignment(group.vector, use.names = TRUE)
    if(length(group)>130){
      #(!) temporary workaround to print partial info for large groups
      print(paste("Warning: group trimmed:",i))
      group.vector <- group.vector[1:130]
      mg <- Biostrings::AAMultipleAlignment(group.vector, use.names = TRUE)
    }
    group.consensus <- msa::msaConsensusSequence(mg, thresh = consensus.thresh)
    write(group.consensus,cpath,append=TRUE)
    num.consensus.pos <- stringr::str_count(group.consensus, "[A-Z]")

    # == == == == == print sequence alignment == == == == ==

    if(length(group) >= min.groupsize & num.consensus.pos >= min.consensus.pos){


      #establish output file names for this group
      pad.num <- formatC(i, width = num.length, flag = "0")
      tname <- paste0(mpath,"msa-", pad.num, ".tex")
      pname <- paste0(mpath,"msa-", pad.num, ".pdf")

      # calculate fig height (inch). top margin 0.75, each line 0.14. key ~1.0?
      msa.height <- 0.75 + 0.14*(length(group.vector) + 1) + 1.3

      #print msa logo
      msa::msaPrettyPrint(mg, output="tex", file = tname,
                     askForOverwrite=FALSE,
                     paperWidth = 12, paperHeight = msa.height,
                     consensusThreshold = rev(consensus.thresh))

      #convert tex to pdf
      if(pdflatex){
        tools::texi2pdf(tname, clean=TRUE)
        file.copy(paste0("msa-",pad.num,".pdf"),pname)
        file.remove(paste0("msa-",pad.num,".pdf"))
        file.remove(tname)


        # #optionally convert pdf to png
        # if(make.png){
        #   iname <- paste0(mpath,"msa-",i,"_1.png")
        #   pdftools::pdf_convert(pname, format = "png", filenames = iname,
        #                         dpi = 300,verbose=FALSE)
        # }


        #cleanup working directory
        file.remove(list.files()[grep("\\.aux|\\.log",list.files())])

      }
    }
  }

  #merge pdfs
  if(pdflatex & pdftk){
    staplr::staple_pdf(input_directory = mpath,
                       output_filepath = paste0(mpath,"msa.pdf"))
  }

}
