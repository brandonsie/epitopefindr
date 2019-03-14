#' groupMSA
#'
#' Generate pdf multiple sequence alignment logos for identified aligning groups.
#'
#' @param groups Table of final peptides and pre-calculated groupings.
#' @param mpath Directory to write sequence alignment files.
#' @param min.groupsize Minimum number of peptides per group to require in order to print.
#' @param pdflatex Logical whether or not to produce PDF LaTeX figures using pdflatex
#' @param pdftk Logical whether or not to use staplr and pdftk to merge individual msa pdfs.
#' @param trim.groups Logical whether or not to apply msaTrim to edges of logos. Not implemented.
#' @param make.png Depreciated. Locial whether or not to convert PDF output to PNG.
#'
#' @export

groupMSA <- function(groups, mpath = "intermediate_files/msa/",
                     min.groupsize = 2, pdflatex = TRUE, pdftk = TRUE,
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

    if(length(group) >= min.groupsize){

      # if(trim.groups){ #optionally trim groups using microseq package
      #   mg <- msa::msaClustalW(group)
      #   mf <- data.frame(Header = names(as.character(mg)),
      #                    Sequence = as.character(mg) %>% as.vector,
      #                    stringsAsFactors = FALSE)
      #
      #   gpath <- paste0(output.dir, "groups/group", i, "_msa.fasta")
      #   microseq::writeFasta(mf, gpath)
      #   mt <-  microseq::msaTrim(microseq::readFasta(gpath), 0, 0)
      #
      #   #convert back to seqinr fasta file
      #   tpath <- paste0(output.dir, "groups/group", i, "_trim.fasta")
      #   seqinr::write.fasta(mt$Sequence %>% as.list, mt$Header, tpath)
      #   group <- Biostrings::readAAStringSet(tpath)
      # }

      # == == == == == normalize peptide name length == == == == ==
      k = 50 #characters from peptide name to use


      #separate basenames and start/end positions
      basenames <- names(group) %>% gsub("\\.[0-9]+\\.[0-9]+$", "", .)
      positions <- stringr::str_extract(names(group), "\\.[0-9]+\\.[0-9]+$")

      #shorten long names
      shortened.basenames <- substr(basenames, 1, k)

      #pad short names
      lengthened.basenames <- stringr::str_pad(shortened.basenames,
                                               width=k+1,side="right",pad=".")

      #reappend start & end positions
      fully.modified.names <- paste0(lengthened.basenames, positions)

      #update group names
      names(group) <- fully.modified.names

      # == == == == == perform sequence alignment == == == == ==
      # cat("\n")
      # mg <- msa::msaClustalW(group) #old way that caused gaps

      #(!) bookmark new way
      group.vector <-as.character(group)
      # group.vector <- groups$Seq[groups$Group == i]
      # names(group.vector) <- groups$ID[groups$Group == i]
      mg <- Biostrings::AAMultipleAlignment(group.vector, use.names = TRUE)

      #(!) temporary workaround to print partial info for large groups
      if(length(group)>130){
        print(paste("Warning: group trimmed:",i))

        # # old way
        # group <- group[1:130]
        # mg <- msa::msaClustalW(group)

        # new way
        group.vector <- group.vector[1:130]
        mg <- Biostrings::AAMultipleAlignment(group.vector, use.names = TRUE)

      }

      write(msa::msaConsensusSequence(mg),cpath,append=TRUE)

      #establish output file names for this group
      pad.num <- formatC(i, width = num.length, flag = "0")
      tname <- paste0(mpath,"msa-", pad.num, ".tex")
      pname <- paste0(mpath,"msa-", pad.num, ".pdf")

      # calculate fig height (inch). top margin 0.75, each line 0.14. key ~1.0?
      msa.height <- 0.75 + 0.14*(length(group.vector) + 1) + 1.3

      #print msa logo
      msa::msaPrettyPrint(mg, output="tex", file = tname,
                     askForOverwrite=FALSE, paperWidth = 12,
                     paperHeight = msa.height)

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
