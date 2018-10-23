# -------- Main Functions ---

epFindExample <- function(){
	epitopeFinder("example","1","any")
}

epitopeFinder <- function(proj.id, e.thresh, g.method = "any", 
												 aln.size = 7, autorun = TRUE){
	# == == == == == Setup/configuration steps. == == == == == 
	if(grepl("\\.fasta",proj.id)){proj.id <- gsub("\\.fasta","",proj.id)} 
	
	print(paste("Running epitope finder on [",proj.id,".fasta] with [e < ",
							e.thresh,"] and grouping method [",g.method,"].",sep=""))
	
	print(paste(format(Sys.time(), "%H:%M:%S"),"Step 1 of 6:",
							"BLASTing input sequences against each other."))
	libcall() #if necessary load relevant R packages
	epSetupDirectory(proj.id, e.thresh, g.method) #prepare output directories
	epSetupPeptides() #cleanup input sequences
	epSetupBLAST() #blast input sequences against each other and prep for analysis
	
	# load some data from global environment
	gl <- c("path","blast.main","blast.id3","blast.id4")
	for(i in gl){assign(i,glGet(i))}
	
	# == == == == == Main script execution. == == == == == 
	if(autorun){
		
		nindex <- blast.main$qID %>% unique %>% length
		print(paste(format(Sys.time(), "%H:%M:%S"),"Step 2 of 6:",
								"Identifying epitopes for",nindex,"peptides."))
		path %<>% pbCycleBLAST(ncycles="max"); fwrite(glGet("blast.main"),blast.id3)
		
		print(paste(format(Sys.time(), "%H:%M:%S"),"Step 3 of 6:",
								"Looping back through peptides in reverse order."))
		path %<>% trimEpitopes(); fwrite(glGet("blast.main"), blast.id4)
		
		print(paste(format(Sys.time(), "%H:%M:%S"),"Step 4 of 6:",
								"Grouping epitope sequences."))
		indexGroups(path, mode = g.method) #grou peptides. mode = "any" | "all"
		
		print(paste(format(Sys.time(), "%H:%M:%S"),"Step 5 of 6:",
								"Generating multiple sequence alignment motifs."))
		groupMSA() #generate multiple sequence alignment motifs
		
		
		print(paste(format(Sys.time(), "%H:%M:%S"),"Step 6 of 6:",
								"Preparing other output files."))
		outputTable() #generate output table
		outputFiles() #copy relevant output files to a new directory
		
		print(paste(format(Sys.time(), "%H:%M:%S"),"epitopeFinder run complete!"))
	}
}

pbCycleBLAST <- function(path, ncycles="max"){
	# == == == == == A. Initialize == == == == ==
	blast.main <- glGet("blast.main")
	path <- glGet("path")
	
	# == == == == == B. Count starting full peptides. == == == == ==
	if(ncycles == "max"){
		n.pep <- blast.main$qID %>% unique %>% length 
		
	} else if(class(ncycles) == "numeric"){
		n.pep <- ncycles
		
	} else {
		stop("Error: pbCycleBLAST: improper ncycles parameter.")
	}
	
	# == == == == == C. Run cycle of epitopeBLAST. == == == == ==
	pb <- epPB(-n.pep,0)
	
	path %<>% cycleBLAST(pb, n.pep)
	glAssign("path", path)
	
	return(path)
	
}

cycleBLAST <- function(path, pb, n, verbose = FALSE){
	
	#Repeatedly call epitopeBLAST until all peptides are converted to epitopes
	if(n > 0){ # n counts down to zero
		setTxtProgressBar(pb, -n) #update progress bar in console
		
		#Optionally print time at each iteration for debugging / speed optimization
		if(verbose == TRUE) {print(paste(Sys.time(), path, "start"))}
		
		#Check that there  are still unprocessed index peptides (name lacks period)
		peptides <- readAAStringSet(path)
		if(names(peptides) %>% grepl("\\.", .) %>% mean < 1){
			
			path <- epitopeBLAST(path, verbose)
			glAssign("path", path)
			path <- glGet("path")
			n <- n - 1 # -sum(!(grepl("\\.", names(epitopesI))))
			
			path <- cycleBLAST(path, pb, n)
		} 
	}
	
	return(path)
}

epitopeBLAST <- function(path, verbose = FALSE){
	# == == == == == A. Load fasta, separate full peptides & epitopes. == =
	peptides <- readAAStringSet(path)
	epitopes <- peptides[grepl("\\.", names(peptides))]
	
	gl <- c("blast.main","output.dir","proj.id")
	for(i in gl){assign(i,glGet(i))}
	
	
	#subset out prior epitopes. ep.genes = gene|tile names of epitopes.
	ep.names <- names(epitopes) %>% strsplit("__") %>% unlist %>% as.character
	ep.genes <- sapply(1:length(ep.names), function(x){
		ep.names[x] %>% strsplit("\\.") %>% unlist %>% `[[`(1)})
	blast.current <- blast.main[!(blast.main$qID %in% ep.genes), ] 
	
	# == == == == == B. Choose index peptide & imsadentify its epitopes. == =
	index <- blast.current %>% chooseIndex()
	ipath <- paste0(output.dir,"epitopes/","indexOrder.txt")
	write(index,ipath,append=TRUE)
	
	blast.index <- rbind(
		blast.main[blast.main$qID==index,-"nAlign"],
		qsSwap(blast.main[blast.main$sID==index,-"nAlign"])) %>% unique
	blast.backup <- blast.main
	
	if(!exists("verbose")){verbose <- FALSE}
	if(verbose == FALSE){
		epIndex <- blast.index %>% indexEpitopes
	} else{
		print(path)
		print(paste(index, ", Blast Rows:", nrow(blast.index)))
		system.time(epIndex <- blast.index %>% indexEpitopes %>% print)
		print(as.character(epIndex[, 1]))
		cat("\n")
	} 
	
	# == == == == == C. Write new epitopes .fasta file. == = 
	#replace index peptides with index epitopes
	blast.remain <- blast.current[blast.current$qID!=index, ] #aln, q != index
	epFrag <- data.frame(ID = names(epitopes), Seq = epitopes %>% as.character)
	epList <- data.frame(ID = blast.remain$qID, Seq = blast.remain$qSeq) %>% 
		rbind(epIndex, epFrag) %>% unique #blast.remain + old eps + new index eps
	epList %<>% mergeFastaDuplicates #remove duplicates, merge nomenclature
	
	#write alignments: index epitopes, remaining peptides, former epitopes
	fpath <- paste0(output.dir, "epitopes/")
	if(!dir.exists(fpath)){dir.create(fpath)}
	if(path == paste0(output.dir, "epitopes/", proj.id, ".fasta")){
		path <- paste0(output.dir, "epitopes/epitopes1.fasta")
	} else{
		gpath <- substr(path, 3, nchar(path)) %>% strsplit("/") %>% unlist
		j <- parse_number(gpath[length(gpath)]) + 1
		path <- paste0(output.dir, "epitopes/epitopes", j, ".fasta")
	}
	writeFastaAA(epList, path)
	return(path)
} #END epitopeBLAST()

indexEpitopes <- function(blast.index){
	# blast.main <- fread(paste0(output.dir, "sjogrens_blast_2_precycle.csv"))
	# == == == == == A. Set up data frames. == =
	pos00 <- blast.index[, c("qStart", "qEnd")] %>% unique %>% data.frame #
	pos00 <- pos00[(pos00$qEnd - pos00$qStart > 5), ]
	pos00 <- pos00[order(pos00$qEnd, pos00$qStart), ]	#primary sort by qEnd!
	rownames(pos00) <- c(1:nrow(pos00))
	
	pos0 <- pos00
	pos1 <- setnames(data.frame(matrix(nrow = 0, ncol = 2)), names(pos00))
	df <- data.frame(findOverlaps(makeIR(pos00), minoverlap = 7))
	prev <- matrix(nrow = nrow(pos00), ncol = 1)
	
	# == == == == == Run overlap identification cycle == =
	for (i in 1:nrow(pos0)) {
		#for each aln, starting with earliest end
		#get list of all overlaps, sorted by query then by start pos of subject
		
		if ((pos0$qEnd[i] - pos0$qStart[i] > 5)) {
			# == == == == == First Epitope == == == == == 
			if (i == 1) {
				#update sdf and psd. this must come after the .pre assignment
				sh0 <- sh.new <- df$subjectHits[df$queryHits == i] %>% sort #aln overlaps
				psd.new <- pos0[sh.new, ] #take corresponding start/end from pos0
				
				#produce minimal epitope
				minep <- c(max(psd.new$qStart), min(psd.new$qEnd))
				pos1 <- setnames(rbind(pos1, minep), names(pos00))
				
				# crop other alns that start before (earliest end-5) to start there
				minst <- min(psd.new$qEnd) - 5
				pos0$qStart[pos0$qStart < minst] <- minst
				pos0[pos0$qEnd - pos0$qStart <= 5, ] <- c(0, 0) #remove from cropped
				
			} else {
				# == == == == == Subsequent Epitopes == == == == == 
				prev[i - 1] <- sh0 %>% list #keep previous data
				sh0 <- sh <- df$subjectHits[df$queryHits == i] %>% sort #new aln overlp
				for(j in sh0[sh0<i]){sh <- sh[!(sh %in% unlist(prev[j]))]} #excl prev
				
				if(length(sh0[sh0<i]) == 0){sh <- c(sh, i) %>% unique %>% sort}
				sh.new <- sh
				
				#take alignments to i that didn't align to previous
				# sh.new <- rem[!(rem %in% c(1:nrow(pos0))[pos0$qEnd == 0])]
				if(length(sh.new)>0){
					sh.new <- c(sh.new, i) %>% unique %>% sort
					psd.new <- pos00[sh.new, ]
					
					# == == == == == produce minimal epitope == == == == == 
					minep <- c(max(psd.new$qStart), min(psd.new$qEnd))
					pos1 <- setnames(rbind(pos1, minep), names(pos00))
					
					# == == == == == crop alns that start before (end-5) == == == == == 
					minst <- min(psd.new$qEnd) - 5
					pos0$qStart[pos0$qStart < minst] <- minst
					
					#remove alignments with cropped length < 7 (incl earliest end aln)
					pos0[pos0$qEnd - pos0$qStart <= 5, ] <- c(0, 0)
					
				}
			}
		}
	}
	
	# Tidy up new overlap table
	posN <- (countOverlaps(makeIR(pos00), makeIR(pos1), minoverlap = 7)) == 0
	w1 <- paste("CAUTION:", 
							"indexEpitopes/posN len > 0 for index", blast.index$qID[1])
	if(length(posN[posN == TRUE])>0) {print(w1)}
	pos1 <- rbind(pos1, pos00[posN, ]) %>% unique
	gpos <- pos1
	
	# == == == == == C. Update blast.main. == =
	#amend blast.main positions to be trimmed to corresponding overlap index ep
	
	blast.main <- glGet("blast.main")
	for(i in 1:nrow(blast.index)){ #for each 
		if(gpos[gpos$qStart == blast.index$qStart[i] & 
						gpos$qEnd == blast.index$qEnd[i],] %>% nrow == 0){ 
			#ignore exact matches
			
			qpos <- c(1:nrow(blast.main))[blast.main$qID == blast.index$qID[i] & 
																			blast.main$sID == blast.index$sID[i] &
																			blast.main$qStart == blast.index$qStart[i] &
																			blast.main$qEnd == blast.index$qEnd[i] &
																			blast.main$sStart == blast.index$sStart[i] &
																			blast.main$sEnd == blast.index$sEnd[i]]
			spos <- c(1:nrow(blast.main))[blast.main$qID == blast.index$sID[i] & 
																			blast.main$sID == blast.index$qID[i] &
																			blast.main$qStart == blast.index$sStart[i] &
																			blast.main$qEnd == blast.index$sEnd[i] &
																			blast.main$sStart == blast.index$qStart[i] &
																			blast.main$sEnd == blast.index$qEnd[i]]
			
			pos <- c(qpos, spos) %>% na.omit %>% as.numeric
			
			
			if(length(pos)>0){
				
				gOverlap <- isOverlapping2(blast.index[i,c("qStart","qEnd")],gpos)
				for(j in gOverlap){ #for each position of gpos that overlaps
					
					dstart <- gpos[j, 1] - blast.index[i, "qStart"] #g-b. positive
					dend <- gpos[j, 2] - blast.index[i, "qEnd"] #g-b. negative
					
					
					#update so that multiple rows are made for each j that aligns
					b.add <- blast.main[pos,]
					for(p in 1:length(pos)){
						b.add[p,"qStart"] %<>% + dstart
						b.add[p,"qEnd"] %<>% + dend
						b.add[p,"sStart"] %<>% + dstart
						b.add[p,"sEnd"] %<>% + dend
					}
					blast.main <- rbind(blast.main, b.add)
					
					
				}
				
				#then remove original qpos spos rows
				blast.main <- blast.main[-pos, ] %>% unique
				
			}
		}
	}
	
	blast.main %<>% removeSmallAln
	blast.main <- blast.main[!duplicated(
		blast.main[,c("qID","sID","qStart","qEnd","sStart","sEnd")]),]
	blast.main %<>% numAlignments()
	glAssign("blast.main", blast.main)
	
	# == == == == == D. Convert to epitope sequences & output. == =
	indexep <- data.frame(ID=character(), Seq=character())
	for(i in 1:nrow(gpos)){
		ID <- paste(blast.index$qID[1] %>% as.character, 
								gpos[i, 1], gpos[i, 2], sep=".")
		Seq <- blast.index$qSeq[1] %>% substr(gpos[i, 1], gpos[i, 2])
		indexep %<>% rbind(setnames(as.data.frame(t(c(ID, Seq))), names(indexep)))
	}
	
	return(indexep)
	
} #end indexEpitopes

trimEpitopes <- function(path, tofilter = FALSE){
	#for debugging
	if(!exists("path")) {path <- gl.path}
	if(!exists("tofilter")) {tofilter <- FALSE}
	
	#load some data from global environment
	
	gl <- c("blast.id3","output.dir")
	for(i in gl){assign(i,glGet(i))}
	
	blast.main <- fread(blast.id3)
	blast.main %<>% prepareBLAST(tofilter)
	glAssign("blast.main", blast.main)
	
	index.order <- fread(paste0(output.dir,"epitopes/indexOrder.txt"),
											 header = FALSE, sep=".")
	
	
	
	#update blast table in reverse order
	pb <- epPB(1,nrow(index.order))
	for(i in nrow(index.order):1){
		setTxtProgressBar(pb,nrow(index.order)-i)
		blast.main <- glGet("blast.main") %>% as.data.table
		index <- index.order[i] %>% as.character
		
		# print(i)
		# print(index)
		
		blast.index <- rbind(
			blast.main[blast.main$qID==index,-"nAlign"],
			qsSwap(blast.main[blast.main$sID==index,-"nAlign"])) %>% unique
		blast.backup <- blast.main
		
		if(nrow(blast.index)>0){
			indexEpitopes(as.data.frame(blast.index))
		}
	}
	
	blast.main <- glGet("blast.main")
	
	#output
	mpath <- paste0(output.dir, "epitopes/final_epitopes.fasta")
	finalep <- blast.main[order(blast.main$qID,blast.main$qStart,blast.main$qEnd),
												c("qID", "qStart", "qEnd", "qSeq")] %>% unique
	finalep$Seq <- sapply(1:nrow(finalep), function(x){
		substr(finalep$qSeq[x], finalep$qStart[x], finalep$qEnd[x])
	})
	finalep$ID <- paste(finalep$qID, finalep$qStart, finalep$qEnd, sep=".")
	writeFastaAA(finalep %>% mergeFastaDuplicates, mpath)
	
	glAssign("path",mpath)
	return(mpath)
}

indexGroups <- function(path, mode="any"){
	
	gl <- c("blast.id4","output.dir")
	for(i in gl){assign(i,glGet(i))}
	
	blast.main <- fread(blast.id4, data.table=FALSE)
	glAssign("blast.main",blast.main)
	
	# == == == == == A. Load final epitope list and blast table(s). == =
	epitopes <- readAAStringSet(path) %>% unmergeFastaDuplicates
	
	#use global blast table
	blast.main <- glGet("blast.main") 
	blast.main <- blast.main[blast.main$qEnd-blast.main$qStart >= 6 & 
													 	blast.main$sEnd-blast.main$sStart >= 6, ]
	
	# == == == == == B. For each epitope, find aligning epitopes. == =
	
	#define data frame from blast.main that represents each epitope's alignments
	b.simp <- setnames(data.frame(matrix(ncol=2, nrow=nrow(blast.main))), c("q", "s"))
	b.simp$q <- sapply(1:nrow(blast.main), function(x){
		paste(blast.main[x, c("qID", "qStart", "qEnd")], collapse=".")})
	b.simp$s <- sapply(1:nrow(blast.main), function(x){
		paste(blast.main[x, c("sID", "sStart", "sEnd")], collapse=".")})
	b.simp <- b.simp[b.simp$q %in% names(epitopes) & 
									 	b.simp$s %in% names(epitopes), ] %>% unique
	
	# == == == == == C. Merge above lists either inclusive or exclusive == =
	
	if(mode == "all"){
		# ---------- "Align to All" = Exclusive Groups (more, tighter) ----------
		#all members of a group must align with all other group members
		
		#define matrix of which fragments align with which fragments via b.sip
		galign <- matrix(0, nrow=length(epitopes), ncol=length(epitopes))
		for(i in 1:nrow(galign)){galign[i, i] <- 1} #set identity to 1
		for(i in 1:nrow(b.simp)){ #populate each cell corresponding to a b.simp row
			r = (1:length(epitopes))[names(epitopes) == b.simp$q[i]]
			c = (1:length(epitopes))[names(epitopes) == b.simp$s[i]]
			galign[r, c] <- 1
		}
		rownames(galign) <- colnames(galign) <- names(epitopes)
		
		#for each peptide, define its minimal group(s)
		aln <- list(Alignments = character())
		addAln <- function(aln, toadd){
			if((aln[[1]] %>% length) == 0){aln[1] <- list(toadd)
			} else{aln[(length(aln)+1)] <- list(toadd)}
			return(aln)
		}
		
		for(i in 1:nrow(galign)){ #
			g1 <- galign[galign[i, ] == 1, galign[i, ] == 1]
			if(mean(g1)!=1){ #if there's misalignment within this group, subdivide
				for(j in 1:nrow(g1)){ #for each peptide in g1 group
					g2 <- j
					for(k in 1:nrow(g1)){ #check all elements, add if it keeps mean = 1
						if(g1[c(g2, k), c(g2, k)]%>%mean == 1){g2 %<>% c(k) %>% unique %>% sort}
					}
					# print(g2)
					aln %<>% addAln(list(rownames(g1)[g2]))
				}
			} else{aln %<>% addAln(list(names(epitopes)[galign[i, ] == 1]))}
		}
		
		#remove duplicate groups and subset groups. order by length
		aln %<>% unique 
		len <- sapply(1:length(aln), function(x){aln[x] %>% unlist %>% length})
		aln <- aln[order(-len)]
		rem <- as.vector(matrix(0, nrow=1, ncol=length(aln)))
		
		for(i in 1:(length(aln)-1)){for(j in (i+1):length(aln)){ #long i, shorter j
			long <- aln[i] %>% unlist
			short <- aln[j] %>% unlist
			if((short %in% long) %>% mean == 1){
				# print(paste(i, j))
				rem[j] <- 1
			}
		}}
		aln <- aln[!rem]
		groups <- aln 
		
	} #end "all"
	
	if(mode == "any"){
		# ---------- "Align to Any" = Inclusive Groups (fewer, looser) ----------
		# all members of a group must align with at least one other group member.
		#merge groups until above condition is satisfied
		
		#define data frame to keep track of each epitope's alignments
		#populate ep.align$Align with lists of aligning epitopes based on b.simp
		ep.align <- setnames(data.frame(matrix(ncol=2, nrow=length(epitopes))), 
												 c("Ep", "Align"))
		ep.align$Ep <- names(epitopes)
		for(i in 1:nrow(ep.align)){
			ep.align$Align[i] <- 
				c(ep.align$Ep[i], b.simp$s[b.simp$q == ep.align$Ep[i]]) %>% list}
		
		k <- 0
		while((ep.align$Align %>% unique %>% unlist %>% 
					 length > length(epitopes))){
			k %<>% +1; if(k > 100) stop("ERROR: stuck in a while loop.")
			
			for(i in 1:nrow(ep.align)){
				eg <- sapply(1:nrow(ep.align), function(x){
					if(ep.align$Ep[i] %in% unlist(ep.align$Align[x])){
						return(1)
					} else return(0)
				}) %>% grep(1, .)
				ep.align$Align[i] <- ep.align[eg, "Align"] %>% unlist %>% 
					unique %>% sort %>% list
			}
		}
		
		#take final unique groups, sorted with longest groups first
		groups <- ep.align$Align %>% unique
		len <- sapply(1:length(groups), function(x) groups[x] %>% unlist %>% length)
		groups <- groups[order(-len)]
	} #end "any"
	
	
	# == == == == == D. Write groups to new fasta files. == =
	gpath <- paste0(output.dir, "groups/")
	if(!dir.exists(gpath)) dir.create(gpath)
	for(i in 1:length(groups)){
		ep <- epitopes[names(epitopes) %in% (groups[i] %>% unlist)]
		ID <- names(ep); Seq <- as.character(ep)
		gdf <- data.frame(ID, Seq)
		writeFastaAA(gdf, paste0(gpath, "group", i, ".fasta"))
	}
	print(paste(length(groups), "groups identified."))
	return(length(groups))
}

groupMSA <- function(trim.groups = FALSE, make.png = FALSE){
	
	#setup output paths and directories
	output.dir <- glGet("output.dir")
	gpath <- paste0(output.dir, "groups/")
	mpath <- paste0(output.dir, "msa/")
	cpath <- paste0(mpath,"consensusSequences.txt")
	if(!dir.exists(mpath)) {dir.create(mpath)}
	if(file.exists(cpath)) {file.remove(cpath)}
	
	#loop through each group
	num <- list.files(gpath) %>% parse_number %>% max
	if(num < 1){break} 
	pb <- epPB(1,num)
	
	for(i in 1:num){
		setTxtProgressBar(pb, i)
		group <- readAAStringSet(paste0(gpath, "group", i, ".fasta"))
		group %<>% unmergeFastaDuplicates
		
		if(length(group)>1){
			
			if(trim.groups){ #optionally trim groups using microseq package
				mg <- msa(group)
				mf <- data.frame(Header = names(as.character(mg)),
												 Sequence = as.character(mg) %>% as.vector,
												 stringsAsFactors = FALSE)
				
				gpath <- paste0(output.dir, "groups/group", i, "_msa.fasta")
				writeFasta(mf, gpath)
				mt <- msaTrim(readFasta(gpath), 0, 0)
				
				#convert back to seqinr fasta file
				tpath <- paste0(output.dir, "groups/group", i, "_trim.fasta")
				write.fasta(mt$Sequence %>% as.list, mt$Header, tpath)
				group <- readAAStringSet(tpath)
			}
			
			# == == == == == normalize peptide name length == == == == ==
			k = 50 #characters from peptide name to use
			norig <- strsplit(names(group),"\\.") %>% unlist %>% matrix(nrow=3) %>% t 
			
			gnames <- substr(norig[,1],1,k)
			
			#shorten long names
			n.long <- str_length(norig[,1]) > k
			gnames[n.long] <- paste0(gnames[n.long],".")
			
			#pad short names
			n.short <- str_length(norig[,1]) <= k
			gnames[n.short] <- str_pad(gnames[n.short],width=k+1,side="right",pad=".")
			
			#append start & end positions
			gpos <- paste(norig[,2],norig[,3],sep = ".")
			names(group) <- paste(gnames,gpos,sep=".")
			
			# == == == == == perform sequence alignment == == == == ==
			cat("\n")
			mg <- msa(group)
			write(msaConsensusSequence(mg),cpath,append=TRUE)
			
			#(!) temporary workaround to print partial info for large groups
			if(length(group)>50){
				group <- group[1:50]
				print(paste("Warning: group trimmed:",i))
				mg <- msa(group)
			}
			
			#establish output file names for this group
			tname <- paste0(mpath,"msa-", i, ".tex")
			pname <- paste0(mpath,"msa-", i, ".pdf")
			
			# calculate fig height (inch). top margin 0.75, each line 0.14. key ~1.0?
			msa.height <- 0.75 + 0.14*(length(group) + 1) + 1.3
			
			#print msa logo
			msaPrettyPrint(mg, output="tex", file = tname,
										 askForOverwrite=FALSE, paperWidth = 8,
										 paperHeight = msa.height)

			#convert tex to pdf
			tools::texi2pdf(tname, clean=TRUE)
			file.copy(paste0("msa-",i,".pdf"),pname)
			file.remove(paste0("msa-",i,".pdf"))
			file.remove(tname)
			
			#optionally convert pdf to png
			if(make.png){
				iname <- paste0(mpath,"msa-",i,"_1.png")
				pdftools::pdf_convert(pname, format = "png", filenames = iname,
										dpi = 300,verbose=FALSE)
			}
			
			#cleanup working directory
			file.remove(list.files()[grep("\\.aux|\\.log",list.files())])
		}
	}
}

outputTable <- function(){
	gl <- c("output.dir","proj.id","path","path0","blast.id4","e.thresh","g.method")
	for(i in gl){assign(i,glGet(i))}
	
	
	#load epitope info
	epitopes <- readAAStringSet(path) %>% unmergeFastaDuplicates
	ep <- data.frame(ID = names(epitopes), Seq = as.character(epitopes))
	original <- readAAStringSet(path0)
	
	## load singletons
	spath <- paste0(output.dir, "epitopes/singletons.csv")
	sing <- fread(spath)
	sf <- data.frame(ID = paste(sing$qID, sing$qStart, sing$qEnd, sep="."), 
									 Seq = sing$qSeq)
	
	#merge epitopes and singletons
	full <- rbind(ep, sf); rownames(full) <- c(1:nrow(full))
	haspipe <- ifelse(length(grep("\\|",full$ID)) > 0,TRUE,FALSE)
	
	#load msa consensus sequences
	cpath <- paste0(output.dir,"msa/consensusSequences.txt")
	cseqs <- gsub("-","",readLines(cpath))
	
	#define and start populating output table
	cnames <- c("gene", "tile", "start", "end", "seq_ep", "seq_tile", 
							"grpnum", "grpname","grpsize","msaconsensus") 
	output <- setnames(data.frame(matrix("NA", nrow=nrow(full), 
																			 ncol=length(cnames))), cnames)
	
	
	if(haspipe){
		output[, c("gene", "tile")] <- strsplit(full$ID, "\\|") %>% 
			unlist %>% as.vector %>% matrix(nrow=2) %>% t
	} else{
		output[, "gene"] <- output[,"tile"] <- full$ID
	}
	
	
	output[, c("tile", "start", "end")] <- strsplit(output$tile, "\\.") %>% 
		unlist %>% as.vector %>% matrix(nrow=3) %>% t
	output$seq_ep <- full$Seq %>% as.character
	
	if(haspipe){
		output$seq_tile <- original[paste(output$gene,output$tile,sep="|")] %>% 
			as.character
	} else {
		output$seq_tile <- original[output$tile] %>% as.character
	}
	
	#load group info
	gpath <- paste0(output.dir, "groups/")
	gmax <- list.files(gpath) %>% parse_number %>% max
	for(i in 1:gmax){
		gi <- readAAStringSet(paste0(gpath, "group", i, ".fasta"))
		gn <- paste(names(gi), collapse="__")
		
		output$grpnum[full$ID %in% names(gi)] %<>% paste(i, sep=", ")
		output$grpsize[full$ID %in% names(gi)] %<>% paste(length(gi), sep=", ")
		output$grpname[full$ID %in% names(gi)] <- gn
	}
	output$grpnum <- gsub("NA, ", "", output$grpnum)
	output$grpsize <- gsub("NA, ", "", output$grpsize)
	
	output$msaconsensus <- suppressWarnings(cseqs[output$grpnum %>% as.numeric])
	
	
	# sort
	output <- output[(order(output$gene,output$tile,
													as.numeric(output$start),
													as.numeric(output$end))),]
	
	opath <- paste0(output.dir, proj.id,"_e-",e.thresh,"_g-",g.method,
									"_outputTable.csv")
	fwrite(output, opath)
	
}

outputFiles <- function(){
	gl <- c("output.dir","proj.id","e.thresh","g.method")
	for(i in gl){assign(i,glGet(i))}
	
	epath <- paste0(output.dir, "epitopes/", proj.id, ".fasta")
	file.copy(epath, paste0(output.dir,"initial_peptides.fasta"))
	
	p1 <- c("Project", "E Threshold", "Group Method")
	p2 <- c(proj.id, e.thresh, g.method)
	pc <- setnames(cbind(p1, p2) %>% data.frame, c("parameter", "value"))
	fwrite(pc, paste0(output.dir, "parameters.txt"))
	
	
	
}

# -------- BLAST Table Manipulators ---

prepareBLAST <- function(blast.thresh, tofilter=TRUE){
	#
	path0 <- glGet("path0")
	fasta <- readAAStringSet(path0)
	blast.merge <- rbind(blast.thresh, qsSwap(blast.thresh)) %>% unique
	blast.tidy <- tidyBLAST(blast.merge, fasta) #update col names, <7aa, gaps
	
	if(!exists("tofilter")){tofilter <- TRUE}
	if(tofilter){
		blast.filter <- filterBLAST(blast.tidy) #remove self-alignments
		rownames(blast.filter) <- c(1:nrow(blast.filter)) #fix row numbering
		return(blast.filter)
	} else{
		rownames(blast.tidy) <- c(1:nrow(blast.tidy))
		return(blast.tidy)
	}
}

tidyBLAST <- function(blast, fasta){
	blast %<>% organizeBLAST() #output table housekeeping (column names, etc.)
	blast %<>% numAlignments() #Add number of alignments per peptide
	blast %<>% addPepSeq(fasta) #add amino acid sequences (tile & align)
	blast %<>% decipherGaps() #split gapped alignments into smaller ungapped
	blast %<>% removeSmallAln() #remove alignmens smaller than 7 aa
	return(blast)
}

filterBLAST <- function(input){ #remove singletons
	output.dir <- glGet("output.dir")
	epath <- paste0(output.dir, "epitopes/")
	spath <- paste0(epath, "singletons.csv")
	
	#remove self alignments
	filter <- input[(input$qID != input$sID) | (grepl("\\.", input$qID)), ] 
	s <- setdiff(input, filter)
	s <- s[!(s$qID %in% filter$qID), ]
	if(!dir.exists(epath)) dir.create(epath)
	if(file.exists(spath)){
		sing <- fread(spath, data.table=FALSE)
		write.csv(rbind(s,sing) %>% unique, spath, row.names = FALSE)
	}
	write.csv(s, spath, row.names = FALSE)
	
	
	return(filter)
}

epSetupDirectory <- function(proj.id, e.thresh, g.method){

	#assign project-specific settings to global environment for easy reference
	glAssign("proj.id", proj.id)
	glAssign("e.thresh", e.thresh)
	glAssign("g.method", g.method)
	
	#create output directories
	output.dir <- paste0("../output/",proj.id,"/",Sys.Date(),"_",
											 format(Sys.time(),"%H%M%S"),"_",
											 proj.id, "_EpitopeFinder/")
	if(!dir.exists(output.dir)){dir.create(output.dir)}
	
	glAssign("output.dir", output.dir)
	if(!dir.exists(output.dir)) {dir.create(output.dir)}
	
	epath <- paste0(output.dir, "epitopes/")
	if(!dir.exists(epath)){dir.create(epath)}
	
	bpath <- paste0(output.dir,"blast/")
	if(!dir.exists(bpath)){dir.create(bpath)}
	
	#prepare directory path references to input .fasta sequences and blast tables
	glAssign("path", paste0("../input/", proj.id, ".fasta")) #starting peptides
	glAssign("blast.id1",paste0(bpath,proj.id,"_blast_1_raw.csv"))
	glAssign("blast.id2",paste0(bpath,proj.id,"_blast_2_precycle.csv"))
	glAssign("blast.id3",paste0(bpath,proj.id,"_blast_3_cycled.csv"))
	glAssign("blast.id4",paste0(bpath,proj.id,"_blast_4_trimmed.csv"))
	
}

epSetupPeptides <- function(){
	#tidy input peptide sequences
	path <- glGet("path")
	output.dir <- glGet("output.dir")
	proj.id <- glGet("proj.id")
	
	fasta <- tidyFasta(path) #remove filler .&* Update path & return AAString
	glAssign("fasta", fasta)
	
	fpath <- paste0(output.dir, "epitopes/", proj.id, ".fasta")
	writeFastaAA(fasta, fpath) 
	glAssign("path", fpath) 
	glAssign("path0", fpath)
}

epSetupBLAST <- function(){
	#run BLASTp on input peptides against each other & pre-process for use
	
	gl <- c("path","blast.id1","blast.id2","e.thresh")
	for(i in gl){assign(i,glGet(i))}
	
	#check for previously written blast table and only re-compute if not present
	if(file.exists(blast.id2)){
		print(paste("Loading existing blast table from", blast.id2))
		blast.main <- fread(blast.id2)
		
	} else {
		if(file.exists(blast.id1)){
			print(paste("Loading existing blast table from", blast.id1))
			blast.path <- fread(blast.id1)
		} else{
			print("Blasting starting sequences against each other...")
			blast.path <- data.table(selfBLASTaa(path)) #blast seq against eachother
			fwrite(blast.path, blast.id1) #write blast to csv 
		}
		
		blast.thresh <- blast.path[blast.path$E < as.numeric(e.thresh), ] 
		blast.main <- prepareBLAST(blast.thresh)
		fwrite(blast.main, blast.id2)
		
	}

	glAssign("blast.main", blast.main)
}

selfBLASTaa <- function(path){ 
	#blastp aa sequences against selves
	makeblastdb(path, dbtype="prot")
	blastdb <- blast(path, type="blastp")
	blastseq <- readAAStringSet(path)
	blastpred <- predict(blastdb, blastseq)
}

queryBLASTaa <- function(query, db){ #input = path to .fasta files
	#blast query sequence against db sequences
	makeblastdb(db, dbtype="prot")
	blastdb <- blast(db, type="blastp")
	blastseq <- readAAStringSet(query)
	blastpred <- predict(blastdb, blastseq)
}

removeSmallAln <- function(blast){
	#delete rows in blast table for which the alignment is fewer than 7aa
	toosmall <- (1:nrow(blast))[(blast$qEnd-blast$qStart)<7]
	if(length(toosmall)>0){blast <- blast[-toosmall, ]}
	return(blast)
}

organizeBLAST <- function(input){
	#rename unwieldy columns from rBLAST output table
	names.orig <- c("QueryID", "SubjectID", "Gap.Openings", 
									"Q.start", "Q.end", "S.start", "S.end")
	names.rep <- c("qID", "sID", "Gaps", 
								 "qStart", "qEnd", "sStart", "sEnd")
	
	#remove unecessary columns 
	for(i in 1:length(names.orig)){names(input)[
		names(input) == names.orig[i]] <- names.rep[i]}
	
	
	names.rem <- c("Perc.Ident", "Alignment.Length", "Mismatches")
	for(i in 1:length(names.rem)){
		if(names.rem[i] %in% names(input)){
			input <- input[, -names.rem[i], with=FALSE]
		}
	}
	# if(ncol(input) == 12)
	# input <- input[, -c(3:5)]
	return(input)
}

decipherGaps <- function(blast){
	#split gapped alignments into multiple smaller ungapped alignments
	
	gpos <- c(1:nrow(blast))[blast$Gaps>0]
	if(length(gpos) == 0){return(blast)}
	
	pb <- epPB(0,length(gpos))	
	pbcount <- 0
	
	
	for(i in gpos){
		pbcount %<>% +1
		setTxtProgressBar(pb, pbcount)
		
		#run pairwise alignment 
		bl.i <- blast[i,]
		qfrag <- with(bl.i,substr(qSeq,qStart,qEnd))
		sfrag <- with(bl.i,substr(sSeq,sStart,sEnd))
		
		msabl <- capture.output(msa(c(qfrag,sfrag) %>% 
															 	as.character %>% AAStringSet))
		msabl <- msabl[(1:2)+length(msabl)-3]
		msagap <- gsub("\\[\\d\\] ","",msabl)
		
		if(nchar(msagap[1]) != nchar(msagap[2])){
			stop("Error: decipherGaps: aligning two sequences with different nchar")}
		
		#make a note of gaps
		g1 <- gregexpr("-", msagap[1])
		g2 <- gregexpr("-", msagap[2])
		g <- c(unlist(g1), unlist(g2))
		
		#identify ungapped fragments
		s <- vector(); e <- vector() #keep track of fragment start/end
		
		if(1 %in% g){gap <- TRUE #check whether starts on a gap
		} else{s <- 1; gap <- FALSE}
		
		for(j in 2:nchar(msagap[1])){ #loop through remaining & make a note of switch
			if(j %in% g){
				#when a new gap starts, mark previous as end
				if(gap == FALSE){e <- c(e, j-1)}
				gap <- TRUE
			} else{
				#when a new non-gap starts, mark as start
				if(gap == TRUE){s <- c(s, j)}
				gap <- FALSE
			}
			if((j == nchar(msagap[1])) & (gap == FALSE)){e <- c(e, j)}
		}
		
		aln <- data.frame(cbind(s, e))
		aln.q <- aln.s <- setnames(data.frame(matrix(nrow=0, ncol=2)), names(aln))
		
		#map fragments to positions in blast alignment table
		for(j in 1:nrow(aln)){
			#need to subtract any gaps that precede an alignment fragment
			g.q <- gregexpr("-", substr(msagap[1], 1, aln[j, 1])) %>% unlist
			g.s <- gregexpr("-", substr(msagap[2], 1, aln[j, 1])) %>% unlist
			
			if(g.q[1] == -1){g.q <- 0} else{g.q %<>% length}
			if(g.s[1] == -1){g.s <- 0} else{g.s %<>% length}
			
			aln.q %<>% rbind(aln[j, ] - g.q)
			aln.s %<>% rbind(aln[j, ] - g.s)
		}
		
		#incrememnt based on alignment start position in bl.i
		aln.q <- aln.q + bl.i$qStart-1
		aln.s <- aln.s + bl.i$sStart-1
		
		#update blast alignment table
		bn <- blast[i, ]
		for(j in 1:nrow(aln)){
			bn[1, c("Gaps", "qStart", "qEnd", "sStart", "sEnd")] <- 
				c(0, aln.q[j, ], aln.s[j, ]) %>% data.frame
			blast %<>% rbind(bn)
		}
	}
	
	blast <- blast[-gpos, ] #remove old rows
	
}

qsSwap <- function(blast.in){
	#to force blast table to be symmetrical, swap subject and query info
	#make a list of which columns refer to query, subject, and either
	blast.in <- as.data.frame(blast.in)
	sNew <- qOld <- grep("^q|^Q",names(blast.in))
	qNew <- sOld <- grep("^s|^S",names(blast.in))
	
	blast.out <- data.frame(matrix(ncol=ncol(blast.in),nrow=nrow(blast.in))) %>%
		setnames(names(blast.in))
	
	for(i in 1:ncol(blast.in)){
		if(i %in% qNew){ #if this column should now be query
			qNewPos <- grep(paste0("^",i,"$"),qNew) #find which variable correspond
			qOldPos <- qOld[qNewPos] #and take that variable from the old position
			blast.out[,i] <- blast.in[,qOldPos]
		} else if(i %in% sNew){
			sNewPos <- grep(paste0("^",i,"$"),sNew)
			sOldPos <- sOld[sNewPos]
			blast.out[,i] <- blast.in[,sOldPos]
		} else{
			blast.out[,i] <- blast.in[,i]
		}
	}
	
	blast.out <- data.table(blast.out)
	
}

numAlignments <- function(input){
	#add/update column to blast table listing # of alignments per query sequence
	input$nAlign <- apply(input, 1, function(x) sum(input$qID == x[1]))
	return(input)
}

addPepSeq <- function(blast, fasta){
	#add/update column in blast table listing peptide names
	if(class(fasta) == "character") fasta <- readAAStringSet(fasta)
	fdf <- data.frame(id = names(fasta), seq = fasta %>% as.character, 
										stringsAsFactors = FALSE)
	
	#add/update column in blast table listing peptide sequence
	for(i in 1:nrow(blast)){
		blast$qSeq[i] <- fdf$seq[fdf$id == blast$qID[i]]
		blast$sSeq[i] <- fdf$seq[fdf$id	 == blast$sID[i]]
	}
	return(blast)
}


# -------- Other Helper Functions ---

libcall <- function(){
	library(data.table) #fread, fwrite
	library(magrittr) #%<>%
	library(seqinr) #write.fasta
	library(msa)
	library(stringr) #str_count
	library(tools) #texi2pdf
	library(pdftools) #pdf_convert
	library(EBImage) #readimage
	library(rBLAST) #blast, predict
	library(readr) #parse_number
	library(microseq) #writeFasta & msaTrim 
	options(stringsAsFactors = FALSE)
}

# glLoad <- function(path){
# 	#load global parameters file from specified path in drive
# 	glpath <- paste0(path,"glparameters.txt")
# 	
# 	#if doesn't exist, start new empty file
# 	assign("gl.parameters", matrix(ncol=2) %>% data.frame %>% 
# 				 	setnames(c("Param","Value")), envir = .GlobalEnv)
# }


glAssign <- function(name, data){
	#Quick wrapper function for assigning data to global environment
	assign(paste0("gl.", name), data, envir = .GlobalEnv)
}

glGet <- function(name){
	get(paste0("gl.", name), envir = .GlobalEnv)
}

makeIR <- function(df){
	#converts a two-column data frame to an IRanges object for overlap calcs
	if(class(df)[1] == "numeric" & length(df) == 2){
		df <- data.frame(df) %>% t
		
	} else{df <- data.frame(df)}
	
	IRanges(df[, 1], df[, 2])
}

isOverlapping2 <- function(pattern, table, almin = 7){
	#checks whether the stard/end pos in pattern overlaps with the
	#start/end pos of each row in table by at least almin
	#returns positions in table that overlap. if no overlap, returns
	
	opos <- findOverlaps(makeIR(pattern), makeIR(table), 
											 minoverlap = almin) %>% subjectHits
	return(opos)
}

writeFastaAA <- function(seq, fpath){
	#takes data frame with $Seq and $ID or AAStringSet and writes fasta file
	
	if(class(seq)[1] == "data.table"){seq <- data.frame(seq)}
	
	if(class(seq) == "data.frame"){
		write.fasta(seq$Seq %>% as.list, seq$ID, fpath, nbchar = 90, as.string=TRUE)
	} else if(class(seq) == "AAStringSet"){
		write.fasta(seq %>% as.character %>% as.vector %>% as.list, 
								names(seq), fpath, as.string=TRUE, nbchar = 90)
	} else stop("ERROR: unrecognized input to writeFastaAA")
}

tidyFasta <- function(input, name){
	#remove filler sequences after "." stop codon and cterminal "*"
	
	if(class(input) == "character"){
		name <- input
		input <- readAAStringSet(input)
	}
	
	fasta.df <- data.frame(ID = (input %>% names %>% as.character), 
												 Seq = (input %>% as.character))
	
	fasta.df$ID %<>% gsub(" ", "_", .)
	fasta.df$ID %<>% gsub("\\.", "_", .)
	fasta.df$ID %<>% gsub("_", "-", .)
	
	fasta.df$Seq <- apply(fasta.df, 1, function(x){
		strsplit(x[2], "\\.")[[1]][1]
	})
	fasta.df$Seq %<>% gsub("\\*", "", .)
	
	write.fasta(as.list(fasta.df$Seq), fasta.df$ID, name)
	output <- readAAStringSet(name)
}

mergeFastaDuplicates <- function(input){
	#if two fasta entries have the same sequence, keep only one copy of the
	# sequence but concatenante the name to be Name1__Name2
	
	if(class(input)[1] == "AAStringSet"){
		input <- data.frame(ID = names(input), Seq = as.character(input))
	}
	
	input <- unique(input)
	if(sum(duplicated(input$Seq))>0){
		dup <- input$Seq[duplicated(input$Seq)] %>% unique
		for(i in 1:length(dup)){
			input$ID[input$Seq == dup[i]] %<>% paste(collapse="__")}
		input <- input[!duplicated(input$Seq), ]
	}
	rownames(input) <- c(1:nrow(input))
	return(input)
}

unmergeFastaDuplicates <- function(input){
	#if a fasta entry has a concatenated name (from a previous call to 
	# mergeFastaDuplicates, then separate the concatenanted name and make
	# a separate entry for each
	np <- input %>% names %>% str_count("__") #number of peptides - 1
	for(k in 1:length(np)){input <- c(input, rep(input[k], np[k]))}
	for(k in 1:length(input)){	#rename duplicated peptides
		names(input)[names(input) == names(input)[k]] <-
			unlist(strsplit(names(input)[k], "__"))}
	return(input)
}

chooseIndex <- function(input, method="min"){ 
	#identify which of a set of remaining peptides from a blast alignment table
	# should serve as the next index peptide. "min" for fewest alignments and
	# "max" for most alignments
	
	#list of remaining full peptide sequences
	fullpep <- input[!grepl("\\.", input$qID), c("qID","nAlign")] %>% unique
	
	if(method == "min"){ #by default, select index peptide with FEWEST alignments
		index <- fullpep$qID[fullpep$nAlign == min(fullpep$nAlign)][1] %>% 
			as.character
	} else if(method == "max"){
		index <- fullpep$qID[fullpep$nAlign == max(fullpep$nAlign)][1] %>% 
			as.character
	} else{stop("Error: chooseIndex: improper selection method.")}
}

dbCleanup <- function(){
	output.dir <- glGet("output.dir")
	
	# remove temporary files made during blast
	epath <- paste0(output.dir, "epitopes/")
	dbf <- list.files(epath)[list.files(epath) %>% grepl("psq|pin|phr", .)]
	dbf <- paste0(output.dir, "epitopes/", dbf)
	file.remove(dbf)
}

epPB <- function(low,high){
	#Brandon's default settings for txtprogressbar
	txtProgressBar(min = low, max = high, initial = low, 
								 char = "=", width = NA, style = 3, file = "")
}
