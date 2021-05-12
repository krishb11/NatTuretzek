# This file contains all the functions related to kallisto you'll need to count the seq data as well as analyse it.
# Only some kinds of visualization is available and you need to decipher your own methods sometimes

# Constructor and Automated function to generate different datatypes   

MyMuse <- function(cdsfile, fastqfiles_dir, replicates) {	
	value <- list(cdsfile=cdsfile, fastqfiles_dir=fastqfiles_dir, num=NULL, tpms=NULL, rpkms=NULL, eff_lengths=NULL, est_counts=NULL, reps=replicates, names=NULL)
	attr(value, "class") <- "MyMuse"
	print("Is it single-end or paired-end? Type 1 for single and 2 for paired and press enter. Same for the rest that will follow")
	value$num <- as.numeric(readline())
	value$names <- RunKallisto(value)
	value$tpms <- GetTPMs(value)
	value$eff_lengths <- Getefflengths(value)
	value$est_counts <- GetFinalCounts(value)
	if(value$num==1) {
		value$rpkms <- GetRPKMs(value)
		value$average_rpkms <- GetAverageRPKMS(value)	
	} else {
		value$fpkms <- GetFPKMs(value)
		value$average_fpkms <- GetAverageFPKMS(value)
	}
	
	return(value)
}

#MyMuse(cdsfile="dmel-all-cds.fastq.gz", fastqfiles_dir="/home/krishna/test_Nat/run_Kallisto/mel", replicates=3)

# Function to get counts from kallisto:
RunKallisto <- function(MuseObject, check=TRUE) {
	kallisto_commands <- c("index", "quant", "bus", "inspect")
	kallisto_file_type <- c("--single","--paired")
	kallisto_mods <- c("-l", "-s")
	bootstrap <- "-b"
	index <- FALSE 
	if(!index) {
		print("Creating indexed cds file.....")
		x <- paste("~/anaconda3/bin/kallisto", kallisto_commands[1], "-i index.idx" , MuseObject$cdsfile, sep=" ")
		if(check) {
			cat(x)
		}
		system(x)
		index <- TRUE
		print("done")
	}
	#Fastq files must be in terms of directory they belong to
	filenames <- list.files(MuseObject$fastqfiles_dir, pattern="*.f*.gz", full.names=TRUE) 
	count <- 1
	print(filenames)
	num <- MuseObject$num	
	if(num==1) {
			print("values of -l and -s respctively")
			l <- as.numeric(readline())
			s <- as.numeric(readline())
			print("How many time you want to do bootstrapping")
			b <- as.numeric(readline())
			names <- lapply(filenames, function(x) {
				x <- unlist(strsplit(x, split="/"))
				x <- tail(x, n=1)
				x <- unlist(strsplit(x, split=".", fixed=TRUE))
				x <- x[1]
				return(x)
			})
			final_output <- lapply(filenames, function(x, kallisto_commands.=kallisto_commands, 
							kallisto_file_type.=kallisto_file_type, kallisto_nums.=MuseObject$kallisto_nums) {
				y <- paste("~/anaconda3/bin/kallisto", kallisto_commands.[2], "-i index.idx -o", names[count], kallisto_file_type.[1], 
					kallisto_mods[1], l, kallisto_mods[2], s, "-b", b, x, sep=" ")
					count <<- count + 1
				return(y)
			})
			lapply(final_output, system)
	} else {
			names <- lapply(filenames, function(x) {
				x <- unlist(strsplit(x, split="/"))
				x <- tail(x, n=1)
				x <- unlist(strsplit(x, split=".", fixed=TRUE))
				x <- x[1]
				x <- unlist(strsplit(x, split="_"))
				x <- x[1]
				return(x)
			})
			names <- as.character(names)
			print("How many time you want to do bootstrapping")
			b <- as.numeric(readline())
			final_output <- lapply(filenames, function(x, kallisto_commands.=kallisto_commands, names.=names, b.=b) {
				y <- paste("~/anaconda3/bin/kallisto", kallisto_commands.[2], "-i index.idx", "-o", names.[count], "-b", b., sep=" ")
				count <<- count + 1
				return(y)
				})
			count <- 1
			final_output <- unlist(final_output)
				for(i in 1:length(final_output)) {
					if(i%%2 == 0) {
						next
					} else {
						z[count] <- paste(final_output[i], filenames[i], filenames[i+1], sep=" ")
						print(z[count])
						system(z[count])
						count <- count + 1
					}
				}
	}
	system('rm -rf index.idx')
	return(names)
}

read_delim <- function(x, sep="\t") {
	return(read.delim(x, sep=sep, header=TRUE))
}

#Go to the directory above all the kallisto directories
GetTPMs <- function(MuseObject) {
	print("acquiring tpms.........")
	filenames <- lapply(MuseObject$names, function(x) {
			x <- list.files(x, pattern="*.tsv", full.names=TRUE)
			return(x)
		})
	files <- lapply(filenames, read_delim)
	files <- lapply(files, as.data.frame)
	#head(files)
	tpms <- lapply(files, function(x) {
			x <- x$tpm
			return(x)
		})
	tpms[[length(tpms)+1]] <- files[[1]][,1]
	names(tpms) <- c(unlist(MuseObject$names), "target_id")
	print("done")
	return(as.data.frame(tpms))	
}

#Go to the directory above all the kallisto directories
Getefflengths <- function(MuseObject) {
	print("acquiring eff_lengths.........")
	filenames <- lapply(MuseObject$names, function(x) {
			x <- list.files(x, pattern="*.tsv", full.names=TRUE)
			return(x)
		})
	files <- lapply(filenames, read_delim)
	files <- lapply(files, as.data.frame)
	eff_length <- lapply(files, function(x) {
			x <- x$eff_length
			return(x)
		})
	eff_length[[length(eff_length)+1]] <- files[[1]][,1]
	names(eff_length) <- c(unlist(MuseObject$names), "target_id")
	print("done")
	return(as.data.frame(eff_length))
}

GetFinalCounts <- function(MuseObject) {
	print("acquiring estimated counts.......")
	filenames <- lapply(MuseObject$names, function(x) {
			x <- list.files(x, pattern="*.tsv", full.names=TRUE)
			return(x)
		})
	files <- lapply(filenames, read_delim)
	files <- lapply(files, as.data.frame)
	est_counts <- lapply(files, function(x) {
			x <- x$est_counts
			return(x)
		})
	est_counts[[length(est_counts)+1]] <- files[[1]][,1]
	names(est_counts) <- c(unlist(MuseObject$names), "target_id")
	print("done")
	return(as.data.frame(est_counts))
}

#The output directory in the below file must be the top level directory of kallisto output directory
GetRPKMs <- function(MuseObject) {
	print("acquiring rpkms......")
	counts <- MuseObject$est_counts[,1:(ncol(MuseObject$est_counts)-1)]
	eff_lengths <- MuseObject$eff_lengths[,1]
	indi_rpkms <- function(counts, .eff_lengths=eff_lengths) {
		for(i in 1:length(counts)) {
			v[i] <- counts[i]/(.eff_lengths[i]/1000 * sum(counts)/1000000)
		}
		return(v)
	}
	rpkms <- apply(counts, indi_rpkms)
	d <- cbind(MuseObject$est_counts[,"target_id"], rpkms)
	colnames(d) <- c("target_id", MuseObject$names)
	return(as.data.frame(d))
	print("done")
}

GetAverageRPKMS <- function(MuseObject) {
	print("Acquiring average rpkms......")
	names <- MuseObject$names
	rpkms <- MuseObject$rpkms
	replicates <- MuseObject$reps
	rpkms$target_id <- NULL
	col <- seq(1, ncol(rpkms)-replicates+1, replicates)
	average_rpkms <- lapply(col, function(x) {
			d <- rpkms[, x:(x+replicates-1)]
			d <- rowSums(d)/replicates
			return(d)
		})
	return(as.data.frame(cbind(as.character(MuseObject$rpkms$target_id), average_rpkms)))
}

GetFPKMs <- function(MuseObject) {
	print("acquiring fpkms......")
	counts <- MuseObject$est_counts[,1:(ncol(MuseObject$est_counts)-1)]
	eff_lengths <- MuseObject$eff_lengths[,2]
	indi_fpkms <- function(counts, eff_lengths=eff_lengths) {
		for(i in 1:length(counts)) {
			v[i] <- counts[i]/(eff_lengths[i]/1000 * sum(counts)/1000000)
		}
		return(v)
	}
	fpkms <- apply(counts, indi_fpkms)
	d <- cbind(MuseObject$est_counts[,1], fpkms)
	colnames(d) <- c("target_id", names)
	print("done")
	return(as.data.frame(d))	
}

GetAverageFPKMS <- function(MuseObject) {
	print("Acquiring average rpkms......")
	names <- MuseObject$names
	rpkms <- MuseObject$rpkms
	replicates <- MuseObject$reps
	rpkms$target_id <- NULL
	col <- seq(1, ncol(rpkms)-replicates+1, replicates)
	average_rpkms <- lapply(col, function(x) {
			d <- rpkms[, x:(x+replicates-1)]
			d <- rowSums(d)/replicates
			return(d)
		})
	return(as.data.frame(cbind(as.character(MuseObject$rpkms$target_id), average_rpkms)))
}

reorder_names <- function(MuseObject, order) {
	return(MuseObject$names[order])
}

heatmap_plotter <- function(MuseObject, which_genes, genes=TRUE, order_for_y, labels, order_for_x, colors) {
	rpkms <- MuseObject$rpkms
	names <- MuseObject$names
	if(!genes) {
		rpkms <- rpkms[rpkms$target_id %in% which_genes[,1], ]
	} else {
		mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = 'ensembl.org')
		t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name"), mart = mart)
		t2g <- dplyr::rename(t2g, transcript=ensembl_transcript_id ,gene_id = ensembl_gene_id, ext_gene = external_gene_name)	
		rpkms$target_id <- vapply(target_id, function(x, data=t2g) {
				if(x %in% t2g$ensembl_transcript_id) {
					var <- t2g[t2g$ensembl_transcript_id %in% x, ]
					var <- var$ensembl_gene_id
					return(var)
				}	
			}, character(1))
		list <- list()
		for(i in 1:length(unique(rpkms$target_id))) {
			list[[i]] <- rpkms[rpkms$target_id %in% unique(rpkms$target_id)[i], ]
			list[[i]]$target_id <- NULL
		}
		indi_adder <- function(lister, v=names) {
			sums <- c()
			for(j in 1:length(v)) {
				dummy <- sum(lister[,v[j]])
				sums <- c(sums, dummy)
			}
			names(sums) <- names
			return(sums)
		}
		list <- lapply(list, function(x) {
				x <- indi_adder(x)
				return(x)
			})
		d <- data.frame()
		for(i in 1:length(list)) {
			d[i,] <- list[[i]]
		}
		rpkms <- cbind(as.character(unique(rpkms$target_id)), d)
		colnames(d) <- c("target_id", names)
	}
	rownames(rpkms) <- rpkms$target_id
	rpkms$target_id <- NULL
	data_for_graph <- rpkms %>%
		rownames_to_column() %>%
  			gather(colname, value, -rowname)
			colnames(dt2) <- c("Gene", "Timepoint", "log2rpkm")
	data_for_graph$Timepoint <- factor(data_for_graph$Timepoint, levels=unique(data_for_graph$Timepoint)[order_for_x])
	
	
}




#6hrs <- mel vs 
#L3 <- all
#24hrs <- suz vs
#48hr <- suz vs

