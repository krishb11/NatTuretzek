
#To acquire the longest isoforms from the set of transcripts.
# 
converter <- function(xstring, data) {
		names <- names(xstring)
		list <- list()
		for(i in 1:length(unique(data$genes))) {
			list[[i]] <- data[data$genes %in% unique(data$genes)[i], ]
			names(list)[i] <- as.character(unique(data$genes)[i])
		}
		final <- list[!(lengths(list)==0)]
	c <- 1
	test <- lapply(final, function(c) {
			dummy <- xstring[names(xstring) %in% as.character(c$ids)]
			d <- as.list(dummy)
			count <- 1
			lul <- c()
			namen <- c()
			for(i in 1:length(d)) {
				if(count==1) {
					lul[count] <- length(d[[i]])
					count <- count +1
				} else {
					if(length(d[[i]]) %in% lul) {
						d[[i]] <- "0"
					} else {
						lul[count] <- length(d[[i]])
						count <- count+1
					}
				}
			}
			d <- DNAStringSet(d[!(lengths(d)==1)])
			c <- c[as.character(c$ids) %in% names(d),]
			return(c)	
		})
	for(i in 1:length(test)) {
	dummy <- xstring[names(xstring) %in% test[[i]]$ids]
	if(length(dummy) > 1) {
			dummy <- dummy[which.max(width(dummy))]
			a <- test[[i]]
			a <- as.data.frame(a)
			a <- a[a$ids %in% names(dummy),]
			test[[i]] <- a
		}	
	}
	lol <- lapply(test, function(a) {
			return(as.character(a$ids))
		})
	lol <- unlist(lol)
	xstring <- xstring[names(xstring) %in% lol]
	print("done")
	return(xstring[unique(names(xstring))])
}

genome_editor <- function(genomefile) {
	genome <- readDNAStringSet(genomefile)
	names <- names(genome)
	names <- strsplit(names, split=" ")
	names <- unlist(lapply(names, function(x) {
			return(x[1])
		}))
	names(genome) <- names
	return(genome)
}


removeDups <- function(xstring_from_geneious) {
	data <- readDNAStringSet(xstring_from_geneious)
	names <- names(data)
	names <- strsplit(names, split=" ")
	names <- unlist(lapply(names, function(a) {
				return(a[1])
		}))
	names(data) <- names
	return(data)
}


nameEditor <- function(xstring, nums, splitters) {
	if(length(splitter)==1 && length(nums)==1) {
		names <- names(xstring)
		names <- strsplit(names, split=splitters)
		names <- unlist(lapply(names, function(x) {
					return(x[nums])
			}))
	} else {
		if(length(nums)==1) {
		names <- names(xstring)
		for(i in 1:length(splitters)) {
				names <- strsplit(names, split=splitters[i])
				names <- unlist(lapply(names, function(x) {
					return(x[nums])
				}))
			}
}


getInteractivePlot <- function(sigfile, counts, net, t2g, time, out, bfilesdir=NULL) {
	t2g <- read.csv(t2g)
	net <- read.csv(net)
	if(is.null(bfilesdir)) {
		sigfile <- read.table(sigfile)
		sigfile <- as.character(sigfile$V1)
		counts <- read.csv(counts)
		rownames(counts) <- counts$target_id
		counts$target_id <- NULL
		if(time %in% "L3") {
			counts <- counts[,1:3]
		} else {
			if(time %in% "6") {
				counts <- counts[,4:6]
			} else {
				if(time %in% "24") {
					counts <- counts[,7:9]
				} else {
					counts <- counts[,10:12]
				}
			}
		}
		counts$final <- (counts[,1] +counts[,2] + counts[,3])/3
		counts$final <- log2(counts$final + 1)
		counts <- as.data.frame(cbind(rownames(counts), counts$final))
		colnames(counts) <- c("target_id", "log2tpms")
		counts <- counts[counts$target_id %in% sigfile, ]
		names <- c()
		for(i in 1:length(counts$target_id)) {
			names[i] <- t2g[t2g$genes %in% counts$target_id[i], "names"]
		}
		counts <- as.data.frame(cbind(names, counts$log2tpms))
		colnames(counts) <- c("name", "log2tpms")
		net <- net[net$gene1 %in% sigfile, ]
		net <- net[net$gene2 %in% sigfile, ]
		go <- data.frame()
		for(i in 1:length(net[,"gene1"])) {
			go[i,"gene1"] <- t2g[t2g$genes %in% net[i,"gene1"], "names"]
		}	
		for(i in 1:length(net[,"gene2"])) {
			go[i,"gene2"] <- t2g[t2g$genes %in% net[i,"gene2"], "names"]
		}	
		head(go)

		g <- graph.data.frame(go)
		counts <- counts[counts$name %in% V(g)$name, ]
		counts$log2tpms <- as.numeric(counts$log2tpms)
		colc <- cut(counts$log2tpms, breaks = c(0, 1.5, 4, 7, 10, 13), include.lowest = TRUE)
		col <- c("black", "blue", "grey", "orange", "red")
		ColourScale <- 'd3.scaleOrdinal()
	            .domain(["black", "blue", "grey", "orange", "red"])
	           .range(["#000000" ,"#0000FF" ,"#BEBEBE" ,"#FFA500" ,"#FF0000"]);'
	    sub <- igraph_to_networkD3(g, group=col[colc])

	    p <- forceNetwork(Links = sub$links, Nodes = sub$nodes, 
	             Source = 'source', Target = 'target', NodeID = 'name', 
	             Group = 'group', colourScale=JS(ColourScale), opacityNoHover = 1, 
	                 fontSize = 12, zoom = T)
	    saveWidget(p, file=out)
	} else {
			listoffiles <- list.files(bfilesdir, pattern=".csv")
			lol <- list()
			lol <- lapply(as.list(listoffiles), function(a) {
					a <- as.data.frame(read.csv(a))
					a <- a[,c("target_id", "qval", "b")]
					return(a)
				})
			data <- data.frame()
			for(i in 1:length(lol)) {
				if(i ==1) {
					data <- lol[[i]]
				} else {
					data <- rbind(data, lol[[i]])
				}
			}
			data <- data[!duplicated(data$target_id),]
			data <- data[data$b > 1.5 | data$b < -1.5,]
			data <- data[data$qval < exp(-5),]
	
			net <- net[net$gene1 %in% data$target_id, ]
			net <- net[net$gene2 %in% data$target_id, ]

			uhoh <- c(net$gene1, net$gene2)
			uhoh <- uhoh[!(duplicated(uhoh))]

			data <- data[data$target_id %in% uhoh,]

			go <- data.frame()
			for(i in 1:length(net[,"gene1"])) {
				go[i,"gene1"] <- t2g[t2g$genes %in% net[i,"gene1"], "names"]
			}	
			for(i in 1:length(net[,"gene2"])) {
				go[i,"gene2"] <- t2g[t2g$genes %in% net[i,"gene2"], "names"]
			}	
			head(go)

			names <- c()
			for(i in 1:length(data$target_id)) {
				names[i] <- t2g[t2g$genes %in% data$target_id[i], "names"]
			}
			final_data <- as.data.frame(cbind(names, data$b))
			colnames(final_data) <- c("name", "log2fold")


			g <- graph.data.frame(go)
			sigfile <- read.csv(sigfile)
			data_for_net <- sigfile[,c("target_id", "qval", "b")]
			data_for_net <- data_for_net[data_for_net$b > 1.5 | data_for_net$b < -1.5,]
			data_for_net <- data_for_net[data_for_net$qval < exp(-5),]
			data_for_net <- data_for_net[,c("target_id", "b")]
			names2 <- c()
			for(i in 1:length(data_for_net$target_id)) {
				names2[i] <- t2g[t2g$genes %in% data_for_net$target_id[i], "names"]
			}
			data_for_net <- as.data.frame(cbind(names2, data_for_net$b))
			colnames(data_for_net) <- c("name", "log2fold")
			data_for_net_final <- data_for_net[data_for_net$name %in% final_data$name, ]
			additional <- final_data[!(final_data$name %in% data_for_net$name), ]
			additional$log2fold <- rep(c(0), each=length(additional$name))
			data_for_net_final <- rbind(data_for_net_final, additional)
			#data_for_net <- data_for_net[data_for_net$name %in% V(g)$name, ]
			data_for_net_final$log2fold <- as.numeric(data_for_net_final$log2fold)
			colc <- cut(data_for_net_final$log2fold, breaks = c(-11, -6, -1.5, 1.5, 6, 11), include.lowest = TRUE)
			col <- c("darkblue", "blue", "grey", "orange", "red")	
			ColourScale <- 'd3.scaleOrdinal()
	            .domain(["blue", "skyblue", "grey", "orange", "red"])
	           .range(["#0000FF" ,"#33FFE9" ,"#BEBEBE" ,"#FFA500" ,"#FF0000"]);'
	        sub <- igraph_to_networkD3(g, group=col[colc])
	        p <- forceNetwork(Links = sub$links, Nodes = sub$nodes, 
	             Source = 'source', Target = 'target', NodeID = 'name', 
	             Group = 'group', colourScale=JS(ColourScale), opacityNoHover = 1, 
	                 fontSize = 12, zoom = T)
	    	saveWidget(p, file=out)
	}
}



# Functions post RNASeq analysis

coolgetter <- function(csvfiles, type) {
	for(i in 1:length(csvfiles)) {
		data <- read.csv(csvfiles[i])
		data <- na.omit(data)
		file <- file(paste(type[i], ".txt", sep=""))
		writeLines(data$target_id, file)
		close(file)
	}
	lit <- list()
	for(i in 1:length(type)) {
		lit[[i]] <- read.table(paste(type[i], ".txt", sep=""))
	}
	for(i in 1:(length(type)-1)) {
		lit[[i]] <- lit[[i]][!(lit[[i]]$target_id %in% lit[[length(type)]]$target_id),]
		file <- file(paste(type[1], "high_interest.txt", sep=""))
		writeLines(data$target_id, file)
		close(file)
	}
}

read_files <- function(pattern) {	
	data <- list.files(pattern=pattern)
	names <- strsplit(data, split=".", fixed=TRUE)
	names <- unlist(lapply(names, function(a) {
			return(a[1])
		}))
	for(i in 1:length(data)) {
		dummy <- read.csv(data[i])
		dummy <- na.omit(dummy)
		print(names[i])
		print(length(dummy$target_id))
		up <- dummy[dummy$b>0,]
		down <- dummy[dummy$b<0,]
		write.csv(dummy, paste(names[i], ".csv", sep=""))
		write.csv(up, paste(names[i], "_up.csv", sep=""))
		write.csv(down, paste(names[i], "_down.csv", sep=""))
	}
}

get_ups_and_downs <- function(pattern) {
	data <- list.files(pattern=pattern)
	names <- strsplit(data, split=".", fixed=TRUE)
	names <- unlist(lapply(names, function(a) {
			return(a[1])
		}))
	for(i in 1:length(data)) {
		dummy <- read.csv(data[i])
		print(names[i])
		print(length(dummy$ext_gene))
		file <- file(paste(names[i], ".txt", sep=""))
		writeLines(dummy$target_id, file)
		close(file)
	}
}

get_ids_only <- function(pattern=".txt") {
	data <- list.files(pattern=pattern)
	names <- strsplit(data, split=".", fixed=TRUE)
	names <- unlist(lapply(names, function(a) {
			return(a[1])
		}))
	for(i in 1:length(data)) {
		dummy <- read.table(data[i])
		dummy <- as.character(dummy$V1)
		assign(names[i], dummy, envir=sys.frame(which=0)) 
	}
	return(names)
} 

getVolcano <- function(sigfile_dir, go) {
	setwd(sigfile_dir)
	files <- list.files(sigfile_dir, pattern=".csv", full.names=FALSE)
	names <- strsplit(files, split=".", fixed=TRUE)
	names <- unlist(lapply(names, function(a) {
			return(a[1])
		}))	
	data <- lapply(files, read.csv)
	go <- read.table(go, header=FALSE)
	names(data) <- names
	for(i in 1:length(data)) {
		data[[i]]$diffexpressed <- "NO"
		data[[i]]$diffexpressed[data[[i]]$b > 0.5 & data[[i]]$qval < 0.05] <- "UP"
		data[[i]]$diffexpressed[data[[i]]$b < -0.5 & data[[i]]$qval < 0.05] <- "DOWN"
		data[[i]]$label <- as.character(data[[i]]$ext_gene)
		for(j in 1:length(data[[i]]$label)) {
			if(data[[i]]$target_id[j] %in% go$V1 & (data[[i]]$diffexpressed[j] %in% "UP" | data[[i]]$diffexpressed[j] %in% "DOWN")) {
				next
			} else {
				data[[i]]$label[j] <- ""
			}
			if(data[[i]]$qval[j]==0) {
				data[[i]]$qval[j] <- sort(unique(data[[i]]$qval))[2]
			}
		}
		pdf(paste(names[i], ".pdf", sep=""))
		p <- ggplot(data=data[[i]][!(data[[i]]$label %in% ""), ], aes(x=b, y=-log10(qval), col=diffexpressed, label=label)) +
        	geom_point(alpha=0.25) + 
        	theme_light() +
        	geom_text_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 25) ,size=4) +
        	scale_color_manual(values=c("blue", "red")) +
        	geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        	geom_hline(yintercept=-log10(0.05), col="red") 

        print(p)
        dev.off()
        print(length(data[[i]][!(data[[i]]$label %in% ""), "label"]))
	}
}