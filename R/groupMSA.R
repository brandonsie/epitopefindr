#' groupMSA
#'
#' Generate pdf multiple sequence alignment logos for identified aligning groups.
#'
#' @param groups Table of final peptides and pre-calculated groupings.
#' @param mpath Direcotry to write sequence alignment files.
#' #' @param pdflatex Logical whether or not to produce PDF LaTeX figures using pdflatex
#' @param pdftk Logical whether or not to use staplr and pdftk to merge individual msa pdfs.
#' @param trim.groups Logical whether or not to apply msaTrim to edges of logos. Not implemented.
#' @param make.png Depreciated. Locial whether or not to convert PDF output to PNG.
#'
#' @export

groupMSA <- function(groups, mpath = "intermediate_files/msa/",
                     pdflatex = TRUE, pdftk = TRUE,
                     trim.groups = FALSE, make.png = FALSE){

  #setup output paths and directories
  cpath <- paste0(mpath,"consensusSequences.txt")
  if(!dir.exists(mpath)) {dir.create(mpath)}
  if(file.exists(cpath)) {file.remove(cpath)}

  #loop through each group
  num <- max(groups$Group)
  if(num < 1){stop("Error: no input groups specified.")}

  for(i in 1:num){
    group <- Biostrings::AAStringSet(groups$Seq[groups$Group == i])
    names(group) <- groups$ID[groups$Group == i]
    group %<>% unmergeFastaDuplicates

    if(length(group)>1){

      # if(trim.groups){ #optionally trim groups using microseq package
      #   mg <- msa::msa(group)
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
      norig <- strsplit(names(group),"\\.") %>% unlist %>% matrix(nrow=3) %>% t

      gnames <- substr(norig[,1],1,k)

      #shorten long names
      n.long <- stringr::str_length(norig[,1]) > k
      gnames[n.long] <- paste0(gnames[n.long],".")

      #pad short names
      n.short <- stringr::str_length(norig[,1]) <= k
      gnames[n.short] <- stringr::str_pad(gnames[n.short],width=k+1,side="right",pad=".")

      #append start & end positions
      gpos <- paste(norig[,2],norig[,3],sep = ".")
      names(group) <- paste(gnames,gpos,sep=".")

      # == == == == == perform sequence alignment == == == == ==
      # cat("\n")
      mg <- msa::msa(group)
      write(msa::msaConsensusSequence(mg),cpath,append=TRUE)

      #(!) temporary workaround to print partial info for large groups
      if(length(group)>130){
        group <- group[1:130]
        print(paste("Warning: group trimmed:",i))
        mg <- msa::msa(group)
      }

      #establish output file names for this group
      tname <- paste0(mpath,"msa-", i, ".tex")
      pname <- paste0(mpath,"msa-", i, ".pdf")

      # calculate fig height (inch). top margin 0.75, each line 0.14. key ~1.0?
      msa.height <- 0.75 + 0.14*(length(group) + 1) + 1.3

      #print msa logo
      msa::msaPrettyPrint(mg, output="tex", file = tname,
                     askForOverwrite=FALSE, paperWidth = 12,
                     paperHeight = msa.height)

      #convert tex to pdf
      if(pdflatex){
        tools::texi2pdf(tname, clean=TRUE)
        file.copy(paste0("msa-",i,".pdf"),pname)
        file.remove(paste0("msa-",i,".pdf"))
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
