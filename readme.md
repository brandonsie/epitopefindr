# EpitopeFinder: Minimal Overlaps from BLAST Alignments 
Version: 1.0.0  
Date: October 23, 2018  
Concept: Ben Larman, Daniel Monaco, Brandon Sie  
Author: Brandon Sie  (contact: brandonsie at gmail)  

# Pipeline Overview: 
The purpose of this tool is to describe the alignments among a set of peptide sequences by reporting the overlaps of each peptide's alignments to other peptides in the set. One can imagine inputting a list of peptides enriched by immunoprecipitation to identify corresponding epitopes. 

This script takes a .fasta file listing peptide sequences of interest and calls BLASTp from within R to identify alignments among these peptides. Each peptide's alignments to other peptides are then simplified to the minimal number of "non overlapping" intervals* of the index peptide that represent all alignments to other peptides reported by BLAST. (*By default, each interval must be at least 7 amino acids long, and two intervals are considered NOT overlapping if they share 6 or fewer amino acids). After the minimal overlaps are identified for each peptide, these overlaps are gathered into aligning groups based on the initial BLAST. For each group, a multiple sequence alignment logo (motif) is generated to represent the collective sequence. Additionally, a spreadsheet is written to list the final trimmed amino acid sequences and some metadata. 

![workflow image] (https://github.com/brandonsie/EpitopeFinder/blob/master/graphics/20181025_EpitopeFindRWorkflow.PNG "EpitopeFind.R Workflow")

# Setup:
1. Install [R (version 3.4.2+)](https://www.r-project.org/).  
2. Install [BLAST+ (version 2.7.1+)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
3. Install the following R packages from CRAN: `tools`, `data.table`, `magrittr`, `seqinr`, `stringr`, `pdftools`, `readr`, `microseq`.  
4. Install the following R packages from Bioconductor: `rBLAST`, `EBImage`, `msa`.  
5. Clone this GitHub repo. Source `EpitopeFinder/epitope_script/`. Add your input fasta file to `EpitopeFinder/input/`. 
6. Call `epitopeFinder(proj.id = [name of your .fasta file])` Output data will be written to `EpitopeFinder/output/`.

----------------------------------------------------------------------
# Changelog
* (Version 1.0.01): Bugfix re fasta file peptide naming (gsub commas etc.)
* 2018-10-23 (Version 1.0.0): Bugfix re: alignment gap removal; UX polishing.
* 2018-10-23 (Version 0.2.1): Output directory bugfixes.
* 2018-10-22 (Version 0.2.0): Github version tracking begins. Vectorized some operations to get rid of for loops.
----------------------------------------------------------------------
# Guide
* all-in-one script can be executed with function `epitopeFinder()`
* to use the provided "example" data, you can run `epitopeFinder(proj.id = "example", e.thresh = 1, g.method = "any")`. This example also has the wrapper function `epFindExample()`

`epitopeFinder` calls a few core functions in order:
1. `epSetupDirectory`,`epSetupPeptides`, and `epSetupBLAST` perform preparatory tidying steps and call blastp from BLAST+ to identify alignments among input peptides.
2. `pbCycleBLAST` cycles through each input peptide to find the overlap of its alignment with other peptides from the input. Nested within a call to `pbCycleBLAST` are calls to `epitopeBLAST`, `indexEpitopes`. 
3. `trimEpitopes` performs a second pass through the identified sequences to tidy alignments.
4. `indexGroups` collects trimmed sequences into aligning groups
5. `groupMSA` creates a multiple sequence alignment motif logo for each group
6. `outputTable` creates a spreadsheet summarizing identified sequences and epitope groups
