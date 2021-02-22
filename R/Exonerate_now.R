library("dplyr")
library("Biostrings")
library("readr")
library("rlist")

genome_editor <- function(WitnessObject) {
	genome <- readDNAStringSet(WitnessObject$genome)
	names <- names(genome)
	names <- strsplit(names, split=" ")
	names <- unlist(lapply(names, function(x) {
			return(x[1])
		}))
	names(genome) <- names
	return(genome)
}

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
	ExonerateNow(value)
	return(value)
}

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

GatherCSV <- function() {
	n <- list.files(path="./csvfiles", pattern="*[0-9].csv")
	return(n)
}

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

RemoveFiles <- function() {
	y <- "rm -rf fastafiles"
	system(y)
	z <- "rm -rf csvfiles"	
	system(z)	
}


