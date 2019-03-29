# epitopefindr <img src="https://brandonsie.github.io/docs/HexStickers_epitopefindr_3.png" align="right" width="120">

![License](https://img.shields.io/badge/licence-GPL--3-blue.svg) 
![Code Size](https://img.shields.io/github/languages/code-size/brandonsie/epitopefindr.svg) 
![Last Commit](https://img.shields.io/github/last-commit/brandonsie/epitopefinder.svg)

The purpose of this package is to describe the [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) alignments among a set of peptide sequences by reporting the overlaps of each peptide's alignments to other peptides in the set. One can imagine inputting a list of peptides enriched by immunoprecipitation (e.g. by [PhIP-seq](https://www.nature.com/articles/s41596-018-0025-6)) to identify corresponding epitopes. 

`epitopefindr` takes a .fasta file listing peptide sequences of interest and calls BLASTp from within R to identify alignments among these peptides. Each peptide's alignments to other peptides are then simplified to the minimal number of "non overlapping" intervals* of the index peptide that represent all alignments to other peptides reported by BLAST. (*By default, each interval must be at least 7 amino acids long, and two intervals are considered NOT overlapping if they share 6 or fewer amino acids). After the minimal overlaps are identified for each peptide, these overlaps are gathered into aligning groups based on the initial BLAST. For each group, a multiple sequence alignment logo (motif) is generated to represent the collective sequence. Additionally, a spreadsheet is written to list the final trimmed amino acid sequences and some metadata. 

![workflow](https://raw.githubusercontent.com/brandonsie/brandonsie.github.io/master/docs/EpitopeFindRWorkflow2c.png)


# Setup:  
1. Install [R (version 3.5+)](https://www.r-project.org/).  
2. Install [BLAST+ (version 2.7.1+)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).  
3. In R console, execute: 
``` r  

# Install Bioconductor packages
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "IRanges", "msa", "S4Vectors"))

# Install Github packages
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("mhahsler/rBLAST")  
devtools::install_github("brandonsie/epitopefindr")

# Load & attach
library(epitopefindr)
```


## Optional (Suggested) Additional Setup : 
_(These are not essential to `epitopefindr`, but are used to generate alignment logo PDFs from the alignment data, which can be valuable visualizations.)_  
1. Install a TeX distribution with `pdflatex`. (e.g. [MiKTeX (version 2.9+)](https://miktex.org)). _(Optional; used to convert multiple sequence alignment TeX files to PDF.)_  
2. Install [pdftk (version 2.02+)](https://www.pdflabs.com/tools/pdftk-server/). _(Optional; used to merge individual PDFs into a single file.)_  
  
----------------------------------------------------------------------  
# Guide

1. Prepare a list of your peptides of interest using one of the following two methods. Either of these can be fed as the first input parameter to `epfind`.  
    * Make a [FASTA file](https://zhanglab.ccmb.med.umich.edu/FASTA/) with peptide names and sequences.
    * Make an `AAStringSet` object of peptides (identifier + sequence) as described in the [Biostrings documentation](http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/Biostrings/html/XStringSet-class.html). 
2. To run a typical `epitopefindr` pipeline, try calling `epfind`:
``` r 
# Basic call
epfind(<path to .fasta>, <path to output dir>)

# Without pdflatex or pdftk
epfind(<path to .fasta>, <path to output dir>, 
        pdflatex = FALSE, pdftk = FALSE)

# More stringent e-value threshold
epfind(<path to .fasta>, <path to output dir>, e.thresh = 0.0001)
``` 

A brief summary of the functions called by `epFind2`:  
  * `pbCycleBLAST` cycles through each input peptide to find the overlap of its alignment with other peptides from the input. Nested within a call to `pbCycleBLAST` are calls to `epitopeBLAST`, `indexEpitopes`. 
  * `trimEpitopes` performs a second pass through the identified sequences to tidy alignments.
  * `indexGroups` collects trimmed sequences into aligning groups
  * `groupMSA` creates a multiple sequence alignment motif logo for each group
  * `outputTable` creates a spreadsheet summarizing identified sequences and epitope groups
  
## For more information, please visit:  
* [epitopefindr's pkgdown website](https://brandonsie.github.io/epitopefindr/)  
* [vignette: epitopefindr Application: Viral Peptides](https://brandonsie.github.io/epitopefindr/articles/epitopefindr_avarda_vignette.html)  
