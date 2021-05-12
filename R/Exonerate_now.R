library("dplyr")
library("Biostrings")
library("readr")
library("rlist")



#' Constructor and automated function to read in genome, query and cds files and generate prerequisite files for the step of annotating a genome using Exonerate.
#' @param blasttype: the type of BLAST you want to perform. All kinds of BLASTing methods are available now
#' @param pathtogenome: Path to the genome file(Note: It must be a fasta file)
#' @param pathtoqueryfile: Path to the query file for BLAST. It could be different from the CDS file and in case both are same, do not provide any value to pathtocdsfile
#' @param species: The species in which we are trying to annotate
#' @param pathtocdsfile: Path to CDS file for exonerate in case it is different from the query file for BLAST, initially set to NULL
#' @returns: MyWitness object with all the input parameters and generates output BLAST files in CSV format and also generates BASH scripts that can be modified and used for the EXONERATE step 
MyWitness <- function(blasttype, pathtogenome, pathtoqueryfile, species, pathtocdsfile=NULL) {
	if(is.null(pathtocdsfile)) {
			value <- list(blast=blasttype, genome=pathtogenome, cdsfile=pathtoqueryfile, species=species, query=pathtoqueryfile)
	
		} else {
			value <- list(blast=blasttype, genome=pathtogenome, cdsfile=pathtocdsfile, species=species, query=pathtoqueryfile)
		}
	system("mkdir csvfiles")
	system("mkdir fastafiles")
	attr(value, "class") <- "MyWitness"
	runBlast(value)
	GetDataForExonerate(value)
	ExonerateNow()
	return(value)
}

#' Function that is a part of the MyWitness automation class used to perform BLAST
#' @param WitnessObject: Object of class MyWitness 
#' @param split: A numeric set to split to the csv files for parallelizing the exonerate step
#' @param additional_args: Currently set to NULL but user can use it to provide additional arguments for BLAST
#' @returns:
runBlast <- function(WitnessObject, split=2000, addtitional_args=NULL) {
	if(WitnessObject$blast %in% "tblastx" || WitnessObject$blast %in% "tblastn") {
		x <- paste(WitnessObject$blast, "-subject", WitnessObject$genome,"-query", WitnessObject$query,"-outfmt 10", "-out csvfiles/index.csv", 
		 "-max_target_seqs 1 -num_threads 10"
		 ,sep=" ")
	} else {
		x <- paste("blastn -task", WitnessObject$blast, "-subject", WitnessObject$genome,"-query", WitnessObject$query,"-outfmt 10", "-out csvfiles/index.csv", 
		 "-max_target_seqs 1 -num_threads 10"
		 ,sep=" ")
	}
	
	print(x)
	system(x)
	print("Index created.....")
	index <- read.csv("csvfiles/index.csv", header=FALSE)
	index <- index[,c(1,2)]
	index <- unique(index)
	colnames(index) <- c("transcript", "contig")
	print(as.data.frame(index))
	index$transcript <- as.character(index$transcript)
	index$contig <- as.character(index$contig)
	print("splitting the files into n csv's.......")
	index <- index %>%
	 	group_by(g=ceiling(row_number()/split)) %>%
	 	do(write_csv(., paste0("csvfiles/", WitnessObject$species, .$g[1], '.csv')))
	print("Done")
}


#' Function that is a part of the MyWitness automation class used to generate all the necessary fasta files
#' @param WitnessObject: Object of class MyWitness 
#' @returns: Nothing but generates fasta files that are necessary for exonerate
GetDataForExonerate <- function(WitnessObject) {
	genome <- readDNAStringSet(WitnessObject$genome)
	genome <- as.list(genome)
	dummy <- DNAStringSet()
	for(i in 1:length(names(genome))) {
		dummy <- DNAStringSet(genome[i])
		writeXStringSet(dummy, paste("fastafiles/",names(genome)[i], ".fasta", sep=""))
	}
	cdsfiles <- readDNAStringSet(WitnessObject$cdsfile)
	cdsfiles <- as.list(cdsfiles)
	dummy <- DNAStringSet()
	for(i in 1:length(names(cdsfiles))) {
		dummy <- DNAStringSet(cdsfiles[i])
		writeXStringSet(dummy, paste("fastafiles/", names(cdsfiles)[i], ".fasta", sep=""))
	}

}

#' Function to acquire all the csv files 
GatherCSV <- function() {
	n <- list.files(path="./csvfiles", pattern="*[0-9].csv")
	return(n)
}

#' Function to generate BASH scripts that can be used to perform EXONERATE 
ExonerateNow <- function() {
	n <- GatherCSV()
	for(i in n) {
		data <- read.csv(paste("csvfiles/", i, sep=""), header=TRUE)
		count <- 1
		filename <- unlist(strsplit(i, split=".", fixed=TRUE))
		filename <- filename[1]
		filename <- paste(filename, ".sh", sep="")
		x <- c()
		for(j in 1:length(data$transcript)) {
			gfile <- paste("fastafiles/", as.character(data$contig[j]), ".fasta", sep="")
			cfile <- paste("fastafiles/", as.character(data$transcript[j]), ".fasta", sep="")
			gfffile <- "suz.gff"
			x[count] <- paste("exonerate --model e2g", "-q", cfile, "-t", gfile, "--softmaskquery yes --softmasktarget yes --bestn 1 --minintron 20 --maxintron 20000  --showalignment false --showtargetgff >>", gfffile, sep=" ")
			count <- count + 1
		} 
		fileConn <- file(filename)
		writeLines(x, fileConn)
		close(fileConn)
	}
}

#' Function to remove the fasta and csv files after the completion of Exonerate
RemoveFiles <- function() {
	y <- "rm -rf fastafiles"
	system(y)
	z <- "rm -rf csvfiles"	
	system(z)	
}


