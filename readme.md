# EpitopeFinder: Minimal Overlaps from BLAST Alignments 
Version: 1.0.06  
Date: December 23, 2018  
Concept: Ben Larman, Brandon Sie, Daniel Monaco  
Author: Brandon Sie  (contact: brandonsie at gmail)  

# Pipeline Overview: 
The purpose of this tool is to describe the alignments among a set of peptide sequences by reporting the overlaps of each peptide's alignments to other peptides in the set. One can imagine inputting a list of peptides enriched by immunoprecipitation to identify corresponding epitopes. 

This script takes a .fasta file listing peptide sequences of interest and calls BLASTp from within R to identify alignments among these peptides. Each peptide's alignments to other peptides are then simplified to the minimal number of "non overlapping" intervals* of the index peptide that represent all alignments to other peptides reported by BLAST. (*By default, each interval must be at least 7 amino acids long, and two intervals are considered NOT overlapping if they share 6 or fewer amino acids). After the minimal overlaps are identified for each peptide, these overlaps are gathered into aligning groups based on the initial BLAST. For each group, a multiple sequence alignment logo (motif) is generated to represent the collective sequence. Additionally, a spreadsheet is written to list the final trimmed amino acid sequences and some metadata. 

![workflow](https://raw.githubusercontent.com/brandonsie/EpitopeFinder/master/man/EpitopeFindRWorkflow.png)

# Setup:
1. Install [R (version 3.4.2+)](https://www.r-project.org/).  
2. Install [BLAST+ (version 2.7.1+)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
3. In R console, execute `devtools::install_github("brandonsie/EpitopeFinder")`, then `library(EpitopeFinder)`
4. Call `epitopeFinder(proj.id = [name of your .fasta file])` Output data will be written to `EpitopeFinder/output/`.

----------------------------------------------------------------------
# Changelog
* 2018-12-23 (Version 1.0.06): Converted project to package (not yet fully documented)  
* 2018-11-08 (Version 1.0.05): Updated some default settings and R environment management
* 2018-11-01 (Version 1.0.04): Added .printSignpost function to tidy epitopeFinder() contents
* 2018-10-31 (Version 1.0.03): Consolidated parameter settings; made some functions internal
* 2018-10-30 (Version 1.0.02): Bugfix re output table seq_tile column
* 2018-10-26 (Version 1.0.01): Bugfix re fasta file peptide naming (gsub commas etc.)
* 2018-10-23 (Version 1.0.00): Bugfix re: alignment gap removal; UX polishing.
* 2018-10-23 (Version 0.2.10): Output directory bugfixes.
* 2018-10-22 (Version 0.2.00): Github version tracking begins. Vectorized some operations to get rid of for loops.
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
