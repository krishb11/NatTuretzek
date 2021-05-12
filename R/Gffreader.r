library("dplyr")
library("readr")
library("Biostrings")

#' Constructor and automated function of class MyDeep to read the gff output after exonerate and extract CDS data
#' @param file: path to gff file.
#' @param species: species that was annotated.
#' @param style: Already set to exonerate, can be modified to use for standard gff files.
#' @param genomefile: Path to the genomefile from which we can perform the extraction of new CDS data.  
#' @returns: A list with the gff data, The exon data and finally the fasta data of all the genes(longest isoforms). 
MyDeep <- function(file, species, style, genomefile) {
	value <- list(pathtofile= file, species=species, data=NULL, 
					edited_gff_object=list(), style=style, fasta=NULL, genomefile=genomefile)
	attr(value, "class") <- "MyDeep"
	value$data <- read_gff(value)
	value$edited_gff_object <- gffeditor(value)
	value$fasta <- gff2fastaprinter(DeepObject=value)
	return(value)
}


#' Function to read gff file
#' @param DeepObject: An object of class MyDeep
#' @param pathto file: Can also be used separately from the automation and in such a case this parameter can be used to provide input gff file. Initially set to NULL.
#' @returns: the data in the gff file as a dataframe
read_gff <- function(DeepObject, pathtofile=NULL) {
	#read gff file content
	if(is.null(pathtofile)) {
		file <- DeepObject$pathtofile
	} else {
		file <- pathtofile
	}
	x <- paste("grep -i exonerate:est2genome", file, "> new.gff", sep=" ")
	system(x)
	x <- paste("grep -v 'source-version' new.gff > final.gff")
	system(x)
	system("rm -rf new.gff")
	file <- "final.gff"
	suppressWarnings(gff.input <- readr::read_delim(file = file, delim = "\t", col_names = FALSE, comment = "#"))
	
	if (ncol(gff.input) > 9)
		stop("The gff file format can not store more than 9 columns!", call. = FALSE)
	
	
	#name standardized columns
	gffNames <- c("seqid", "source", "type", "start", "end", "score",
				  "strand", "phase", "attribute")
	
	names(gff.input)[seq(ncol(gff.input))] <- gffNames[seq(ncol(gff.input))]
	
	return(gff.input)
}	

#data <- read_gff("~/Downloads/inshallah/lol3.gff")

# Function to edit the gff file with MyDeep Object
# @param DeepObejct: An object of class MyDeep.
gffeditor <- function(DeepObject) {	
	
	#Filetering out the data we need. Like gene, cds, intron etc
	style <- DeepObject$style
	#newdata <- DeepObject$data %>%
	#            filter(type==which_one)

	#Modifying the names accordingly
	if(style=="exonerate") {
	index <- DeepObject$data$type %in% "gene"
	count <- 1
	d <- c()
	for(i in 1:length(index)) {
		if(index[i]) {
			d[count] <- i
			count <- count + 1
		}
		else {
			next
		}
	}
	list <- list()
	count <- 1
	for(i in 1:length(d)) { 
		if(d[i] == tail(d, n=1)) {
			list[[count]] <- DeepObject$data[d[i]:length(DeepObject$data$type), ]
		}
		else {
			list[[count]] <- DeepObject$data[d[i]:d[i+1], ]
			count <- count + 1
		}
	}
	newdata <- DeepObject$data %>%
	            filter(type=="gene")
	names <- strsplit(newdata$attribute, split=";")
	names <- lapply(names, function(x) {
		x <- x[2]
		return(x)
		})
	names <- unlist(names)
	names <- strsplit(names, split=" ")
	names <- lapply(names, function(x) {
		  x <- x[3]
		  return(x)
		})
	names <- unlist(names)
	
	names(list) <- names 
	#print(list)
	} else {
		stop("code for Standard gff is in  process", .call=FALSE)
	}
	#write_delim(newdata, path=paste(DeepObject$speices, "_edited", ".gff", sep=""), delim="\t", col_names=TRUE)
	list <- lapply(list, function(x) {
			y <- x %>%
				filter(type=="exon")
			return(y)
		})
	#print(list)
	return(list)
}

#' Function to print fasta file of the newly acquired CDS co-ordinates
#' @param DeepObject: An object of class MyDeep.
#' @param edited_gff_onject: Currently set to NULL, but in case you have an edited GFF data, it can be directly used instead.
gff2fastaprinter <- function(DeepObject, edited_gff_object=NULL) {
	fasta <- DNAStringSet()
	if(!(is.null(DeepObject))) {
		edited_gff_object <- DeepObject$edited_gff_object
	} else {
		edited_gff_object <- edited_gff_object
	}
	#genomefile <- DeepObject$genomefile
	genome <- readDNAStringSet(DeepObject$genomefile)
	names <- names(genome)
	names <- strsplit(names, split=" ")
	names <- unlist(lapply(names, function(x) {
			return(x[1])
		}))
	names(genome) <- names
	#print(names)
	#head(DeepObject$edited_gff_object)
	print("starting sequence extraction")
	fastac <- DNAStringSet()
	newdata <- lapply(edited_gff_object, function(x) {
			for(i in 1:length(x$type)) {
				if(x$start[i] < width(genome[names(genome) %in% x$seqid[1]])) {
					fastac[i] <- subseq(genome[names(genome) %in% x$seqid[1]], start=x$start[i], end=x$end[i])	
				}
			}
			#print(fastac)
			return(as.list(fastac))
		})
	#return(newdata)
	print("done with extraction now merginig....")
	set <- newdata
	dat <- c()
	count <- 1
	for(i in 1:length(set)) {
		lol <- as.character(DNAStringSet(set[[i]][[1]]))
		if(length(set[[i]]) > 1) {
				for(j in 2:length(set[[i]])) {
				dummy <- as.character(DNAStringSet(set[[i]][[j]]))
				lol <-  paste(lol, dummy, sep="")
			}
		}
		dat[count] <- lol
		count <- count + 1
	}
	print("done merginig....")
	dat <- DNAStringSet(dat)
	names(dat) <- names(DeepObject$edited_gff_object)
	writeXStringSet(dat, paste(DeepObject$species, "transcripts", "fasta", sep="."))
	return(dat)
}








#####lol#######


# list <- list()
# count <- 1
# for(i in 1:length(d)) { 
# 	if(d[i] == tail(d, n=1)) {
# 		list[[count]] <- da[d[i]:length(da$type), ]
# 	}
# 	else {
# 		list[[count]] <- da[d[i]:(d[i+1]-1), ]
# 		count <- count + 1
# 	}
# }

# names <- lapply(data_mrna, function(x) {
# 		x <- strsplit(as.character(x[9]), split=";")
# 		x <- unlist(lapply(x, function(x) {
# 				return(tail(x, n=2))
# 			}))
# 		x <- strsplit(x, split="=", fixed=TRUE)
# 		x <- unlist(lapply(x, function(x) {
# 				return(x[2])
# 			}))
# 		return(x)
# 	})
# count <- 1
# for(i in 1:length(list)) {
# 	if(!("intron" %in% list[[i]]$type)) {
# 		count <- count + 1
# 	}
# }