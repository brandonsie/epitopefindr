#' epSetupDirectory
#' Setup epitopefinder directory paths
#' @param proj.id Name of peptide cohort.
#' @param e.thresh Maximum allowed e-value for BLAST alignment to be considered.
#' @param g.method Grouping method "any" or "all". See ?indexGroups
#' @export

epSetupDirectory <- function(proj.id, e.thresh, g.method){

  #assign project-specific settings to global environment for easy reference
  glParamAssign("proj.id", proj.id)
  glParamAssign("e.thresh", e.thresh)
  glParamAssign("g.method", g.method)

  #create output directories
  project.dir <- paste0("../output/",proj.id,"/")
  output.dir <- paste0(project.dir, format(Sys.Date(),"%Y%m%d"),
                       "_",format(Sys.time(),"%I%M%S"),
                       (format(Sys.time(),"%p") %>% substr(1,1)),
                       "_",proj.id, "_EpF/")

  if(!dir.exists(project.dir)){dir.create(project.dir)}
  if(!dir.exists(output.dir)){dir.create(output.dir)}
  glParamAssign("output.dir", output.dir)

  epath <- paste0(output.dir, "epitopes/")
  if(!dir.exists(epath)){dir.create(epath)}

  bpath <- paste0(output.dir,"blast/")
  if(!dir.exists(bpath)){dir.create(bpath)}

  #prepare directory path references to input .fasta sequences and blast tables
  glParamAssign("path", paste0("../input/", proj.id, ".fasta"))
  glParamAssign("blast.id1",paste0(bpath,proj.id,"_blast_1_raw.csv"))
  glParamAssign("blast.id2",paste0(bpath,proj.id,"_blast_2_precycle.csv"))
  glParamAssign("blast.id3",paste0(bpath,proj.id,"_blast_3_cycled.csv"))
  glParamAssign("blast.id4",paste0(bpath,proj.id,"_blast_4_trimmed.csv"))

}
