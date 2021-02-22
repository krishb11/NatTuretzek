

MyBeret <- function(pathtocsv, species, conditions, MuseObject=NULL) {
		if(is.null(MuseObject)) {
			value <- list(path=pathtocsv, species=species, conditions=conditions, filteredData=NULL, clusterdata=NULL, genes2keep=NULL, clusters=NULL, data=NULL)
		} else {
			value <- list(data=MuseObject, species=species, conditions=conditions, filteredData=NULL, clusterdata=NULL, genes2keep=NULL, clusters=NULL)
		}
		attr(value, "class") <- "MyBeret"
		value$filteredData <- FilterNow(value)
		value$clusterdata <- ClusterNow(value)
		value$genes2keep <- identifier(value)
		value$clusters <- clusterer(list=datforplot(value))
		write.csv(value$clusters, paste(value$species, "_clusters.csv", sep=""))
	
}

FilterNow <- function(BeretObject) {
	if(is.null(BeretObject$data)) {
		counts <- read.csv(BeretObject$path, header=TRUE)
	} else {
		counts <- MuseObject$est_counts
	}
	rownames(counts) <- counts$target_id
	counts <- counts[,2:ncol(counts)]
	counts <- round(counts)
	dataframe <- data.frame(row.names=colnames(counts), conditions=BeretObject$conditions)
	filter.time <- HTSFilter(counts, conditions, norm="DESeq")
	final_data <- filter.time$filteredData
	final_data <- final_data+1
	final_data <- cbind(as.character(rownames(counts), final_data))
	rownames(final_data) <- NULL
	final_data <- as.data.frame(final_data)
	colnames(final_data)[1] <- "genes"
	write.csv(final_data, paste(species, "_filtered_data", ".csv"))
	return(final_data)
}

ClusterNow <- function(BeretObject) {
	data <- BeretObject$filteredData
	genes <- data$genes
	data$genes <- NULL
	data <- as.matrix(data)
	PMM <- PoisMixClusWrapper(y=data, gmin=1, gmax=25, conds=BeretObject$conditions, norm="DESeq")
	list <- list()
	list$BIC <- data_collector(PMM$BIC.results, BeretObject)
	list$ICL <- data_collector(PMM$ICL.results, BeretObject)
	list$Djump <- data_collector(PMM$Djump.results, BeretObject)
	list$DDSE <- data_collector(PMM$DDSE.results, BeretObject)
	return(list)
}

data_collector <- function(model, BeretObject) {
	result_type <- BeretObject$species
	genes <- as.character(BeretObject$filteredData$genes)
	pdf(file=paste(result_type, "_plots",".pdf", sep=""))
	#plot(model, graphs="lambda")
	plot(model, graphs="map")
	plot(model, graphs="map.bycluster")
	dev.off()

	labels <- cbind(as.character(genes), model$labels)
	colnames(labels) <- c("genes", "labels")
	time <- model$lambda[,1]
	lambda <- cbind(time, model$lambda[,2:ncol(model$lambda)])
	probaPost <- cbind(as.character(genes), model$probaPost)
	write.csv(labels, file=paste(result_type, "_labels",".csv", sep=""))
	write.csv(lambda, file=paste(result_type, "_lambda",".csv", sep=""))
	write.csv(probaPost, file=paste(result_type, "_probaPost",".csv", sep=""))
	list <- list(labels=labels, lambdas=lambda, probaPosts=probaPost)
	return(list)

}

filereader <- function(species) {
	probab_bic <- read.csv(paste(species, "_BIC", "_probaPost.csv", sep=""))
	probab_icl <- read.csv(paste(species, "_ICL", "_probaPost.csv", sep=""))
	probab_DDSE <- read.csv(paste(species, "_DDSE", "_probaPost.csv", sep=""))
	probab_Djump <- read.csv(paste(species, "_Djump", "_probaPost.csv", sep=""))

	labels_bic <- read.csv(paste(species, "_BIC", "_labels.csv", sep=""))
	labels_icl <- read.csv(paste(species, "_ICL", "_labels.csv", sep=""))
	labels_DDSE <- read.csv(paste(species, "_DDSE", "_labels.csv", sep=""))
	labels_Djump <- read.csv(paste(species, "_Djump", "_labels.csv", sep=""))

	list <- list(probab_bic=probab_bic, probab_icl=probab_icl, probab_DDSE=probab_DDSE, 
		probab_Djump=probab_Djump, labels_bic=labels_bic, labels_icl=labels_icl, labels_DDSE=labels_DDSE, labels_Djump=labels_Djump)
	return(list)
}

identifier <- function(BeretObject, file=FALSE) {
	if(!(file)) {
			probab_list <- list(probab_bic=BeretObject$list$BIC$probaPosts, 
						    probab_icl=BeretObject$list$ICL$probaPosts,
					 	    probab_DDSE=BeretObject$list$DDSE$probaPosts,
					 	    probab_Djump=BeretObject$list$Djump$probaPosts)
	} else {
		probab_list <- filereader(species=BeretObject$species)
		probab_list <- probab_list[c("probab_bic", "probab_icl", "probab_DDSE", "probab_Djump")]
	}
	genes2keep <- lapply(probab_list, getter)
	return(genes2keep)
}

getter <- function(dataframe) {
	return(apply(dataframe, 1, function(x) {
			for(i in 2:length(x))
			{
				if(x[i] > 0.99)
				{
					return(x[1])
				}
			}
		}))
}

datforplot <- function(BeretObject) {
	list <- list(labels_bic=BeretObject$list$BIC$labels,
				 labels_icl=BeretObject$list$ICL$labels,
				 labels_DDSE=BeretObject$list$DDSE$labels,
				 labels_Djump=BeretObject$list$Djump$labels)
	genes2keep <- BeretObject$genes2keep
	for(i in 1:4) {
		list[i][[1]] <- list[c(i)][[1]][as.character(list[c(i)][[1]]$genes) %in% as.character(genes2keep[c(i)][[1]]), ] 
		colnames(list[c(i)][[1]]) <- c("genes", "cluster")
	}
	return(list)
}

clusterer <- function(list) {
	for(i in 1:4) {
		clusters <- lapply(list, getter2)
	}
	return(clusters)

}

getter2 <- function(element) {
	newdat <- list()
	count <- 1
	for(j in 1:length(unique(element$cluster))) {
		newdat[[paste("cluster", j)]] <- as.character(element$genes[element$cluster==j])
	}
	return(newdat)
}

lambdaplotter <- function(BeretObject, unique_conditions, palette=NULL, file=FALSE) {
	if(file) {
		lambda_bic <- read.csv(paste(species, "_BIC", "_lambda.csv", sep=""))
		lambda_icl <- read.csv(paste(species, "_ICL", "_lambda.csv", sep=""))
		lambda_DDSE <- read.csv(paste(species, "_DDSE", "_lambda.csv", sep=""))
		lambda_Djump <- read.csv(paste(species, "_Djump", "_lambda.csv", sep=""))
	} else {
		lambda_bic <- BeretObject$list$BIC$lambdas
		lambda_icl <- BeretObject$list$ICL$lambdas
		lambda_DDSE <- BeretObject$list$DDSE$lambdas
		lambda_Djump <- BeretObject$list$Djump$lambdas
	}

	lambda_bic <- lambda_bic[lambda_bic$time[unique_conditions], ]
	lambda_icl <- lambda_icl[lambda_icl$time[unique_conditions], ]
	lambda_DDSE <- lambda_DDSE[lambda_DDSE$time[unique_conditions], ]
	lambda_Djump <- lambda_Djump[lambda_Djump$time[unique_conditions], ]

	pdf(paste("lamdaplots_bic_", species,".pdf"))
	par(mfrow=c(5,5))
	yrange <- c(0,2,4,6)
	colors <- rainbow(ncol(lambda_bic)-1)
	for(j in 2:ncol(lambda_bic)) {
		plot(x=c(1:4), y=yrange, main=paste("Cluster.", j, "no.;",length(clusters$labels_bic[[j]]), sep=" "), type="n", xlab="timepoint", ylab="lambda", xaxt="n",
			cex.lab=0.5, cex.axis=0.5, cex.main=0.5)
		axis(1, at=1:4, labels=unique_conditions, cex.axis=0.5, las=2)
		lines(x=1:4,y=lambda_bic[,j], pch=16, col=colors[j])
	}
	dev.off()

	pdf(paste("lamdaplots_icl_", species,".pdf"))
	par(mfrow=c(5,5))
	yrange <- c(0,2,4,6)
	colors <- rainbow(ncol(lambda_icl))
	for(j in 1:ncol(lambda_icl)) {
		plot(x=c(1:4), y=yrange, main=paste("Cluster.", j, "no.:",length(clusters$labels_icl[[j]]), sep=" "), type="n", xlab="timepoint", ylab="lambda", xaxt="n",
			cex.lab=0.5, cex.axis=0.5, cex.main=0.5)
		axis(1, at=1:4, labels=unique_conditions, cex.axis=0.5, las=2)
		lines(x=1:4,y=lambda_icl[,j], pch=16, col=colors[j])
	}
	dev.off()

	pdf(paste("lamdaplots_DDSE_", species,".pdf"))
	par(mfrow=c(5,5))
	yrange <- c(0,2,4,6)
	colors <- rainbow(ncol(lambda_DDSE))
	for(j in 1:ncol(lambda_DDSE)) {
		plot(x=c(1:4), y=yrange, main=paste("Cluster.", j, "no.:",length(clusters$labels_DDSE[[j]]), sep=" "), type="n", xlab="timepoint", ylab="lambda", xaxt="n",
			cex.lab=0.5, cex.axis=0.5, cex.main=0.5)
		axis(1, at=1:4, labels=c("24hrs", "48hrs", "6hrs", "L3"), cex.axis=0.5, las=2)
		lines(x=1:4,y=lambda_DDSE[,j], pch=16, col=colors[j])
	}
	dev.off()

	pdf(paste("lamdaplots_Djump_", species,".pdf"))
	par(mfrow=c(5,5))
	yrange <- c(0,2,4,6)
	colors <- rainbow(ncol(lambda_Djump))
	for(j in 1:ncol(lambda_Djump)) {
		plot(x=c(1:4), y=yrange, main=paste("Cluster.", j, "no.:",length(clusters$labels_Djump[[j]]), sep=" "), type="n", xlab="timepoint", ylab="lambda", xaxt="n",
			cex.lab=0.5, cex.axis=0.5, cex.main=0.5)
		axis(1, at=1:4, labels=c("24hrs", "48hrs", "6hrs", "L3"), cex.axis=0.5, las=2)
		lines(x=1:4,y=lambda_Djump[,j], pch=16, col=colors[j])
	}
	dev.off()
}

dataframer <- function(data) {
	data$L3 <- (data[,2]+data[,3]+data[,4])/3
	data$ovi6 <- (data[,5]+data[,6]+data[,7])/3
	data$ovi24 <- (data[,8]+data[,9]+data[,10])/3
	data$ovi48 <- (data[,11]+data[,12]+data[,13])/3
	data[,2:13] <- NULL
	return(data)
}

expressionplotter <- function(BeretObject, file=TRUE) {
	if(file) {
		lambda_bic <- read.csv(paste(species, "_BIC", "_lambda.csv", sep=""))
		lambda_bic$time <- NULL
		lambda_icl <- read.csv(paste(species, "_ICL", "_lambda.csv", sep=""))
		lambda_icl$time <- NULL
		lambda_DDSE <- read.csv(paste(species, "_DDSE", "_lambda.csv", sep=""))
		lambda_DDSE$time <- NULL
		lambda_Djump <- read.csv(paste(species, "_Djump", "_lambda.csv", sep=""))
		rownames(lambda_Djump) <- lambda_Djump$time
		lambda_Djump$time <- NULL 
		lambda_Djump <- lambda_Djump[c("L3","6hrs", "24hrs", "48hrs"), ]
		lambda_Djump <- apply(lambda_Djump, 2, function(x) {x <- log2(x+1) })	
	} else {
		lambda_bic <- BeretObject$list$BIC$lambdas

		lambda_icl <- BeretObject$list$ICL$lambdas
		lambda_DDSE <- BeretObject$list$DDSE$lambdas
		lambda_Djump <- BeretObject$list$Djump$lambdas
		lambda_bic$time <- NULL
		lambda_icl$time <- NULL
		lambda_DDSE$time <- NULL
		rownames(lambda_Djump) <- lambda_Djump$time
		lambda_Djump$time <- NULL
		lambda_Djump <- apply(lambda_Djump, 2, function(x) {x <- log2(x+1) })
	}
	clusters <- BeretObject$clusters
	species <- BeretObject$species
	pdf(paste("Expression_plot _BIC_", species, ".pdf", sep=""))
	par(mfrow=c(3,3))
	xrange <- range(1:4)
	yrange <- range(1:200000)
	colors <- rainbow(ncol(lambda_bic))
	for(i in 1:ncol(lambda_bic)) {
		data_new <- data[data$genes %in% clusters$labels_bic[[i]], ]
		data_new[,2:5] <- data_new[,2:5] +1
		heading = paste("Cluster",i, " : ", length(clusters$labels_bic[[i]]), " genes")
  		plot(xrange, yrange, type="n", xaxt = "n", xlab="timepoint", ylab="counts", bty="L", main=heading,
  		cex.lab=1, cex.axis=0.5, cex.main=1, log="y") #first draw an empty plot
  		axis(1, at=1:4, labels=c("L3","6hrs", "24hrs","48hrs"), cex.axis=0.5, las=2)

  		for(j in 1:nrow(data_new)) {
  			lines(x=1:4, data_new[j,2:5], type="b", lwd=0.5, col=colors[i])
  			
  		}

	}
	dev.off()
	pdf(paste("Expression_plot _ICL_", species, ".pdf", sep=""))
	par(mfrow=c(3,3))
	xrange <- range(1:4)
	yrange <- range(1:200000)
	for(i in 1:ncol(lambda_icl)) {
		data_new <- data[data$genes %in% clusters$labels_icl[[i]], ]
		data_new[,2:5] <- data_new[,2:5] +1
		heading = paste("Cluster",i, " : ", length(clusters$labels_icl[[i]]), " genes")
  		plot(xrange, yrange, type="n", xaxt = "n", xlab="timepoint", ylab="counts", bty="L", main=heading,
  		cex.lab=1, cex.axis=0.5, cex.main=1, log="y") #first draw an empty plot
  		axis(1, at=1:4, labels=c("L3","6hrs", "24hrs","48hrs"), cex.axis=0.5, las=2)

  		for(j in 1:nrow(data_new)) {
  			lines(x=1:4, data_new[j,2:5], type="b", lwd=0.5, col=colors[i])
  		}
	}
	dev.off()
	pdf(paste("Expression_plot _DDSE_", species, ".pdf", sep=""))
	par(mfrow=c(5,5))
	xrange <- range(1:4)
	yrange <- range(1:200000)
	for(i in 1:ncol(lambda_DDSE)){
		data_new <- data[data$genes %in% clusters$labels_DDSE[[i]], ]
		data_new[,2:5] <- data_new[,2:5] +1
		heading = paste("Cluster",i, " : ", length(clusters$labels_DDSE[[i]]), " genes")
	  	plot(xrange, yrange, type="n", xaxt = "n", xlab="timepoint", ylab="counts", bty="L", main=heading,
	  	cex.lab=1, cex.axis=0.5, cex.main=1, log="y") #first draw an empty plot
	  	axis(1, at=1:4, labels=c("L3","6hrs", "24hrs","48hrs"), cex.axis=0.5, las=2)
	
	  	for(j in 1:nrow(data_new)) {
	  		lines(x=1:4, data_new[j,2:5], type="b", lwd=0.5, col=colors[i])
	  	}
	}
	dev.off()
	pdf(paste("Expression_plot _Djump_dualy_", species, ".pdf", sep=""))
	par(mfrow=c(3,3))
	xrange <- range(1:4)
	yrange <- range(1:200000)
	zrange <- range(0.05:2.85)
	for(i in 1:ncol(lambda_Djump)) {
		data_new <- data[data$genes %in% clusters$labels_Djump[[i]], ]
		data_new[,2:5] <- data_new[,2:5] +1
		heading = paste("Cluster",i, " : ", length(clusters$labels_Djump[[i]]), " genes")
		head(data_new)
		par(mar = c(5, 4, 4, 4) + 0.3)
  		plot(xrange, yrange, type="n", xaxt = "n", xlab="timepoint", ylab="counts", bty="L", main=heading,
  		cex.lab=1, cex.axis=0.5, cex.main=1, log="y") #first draw an empty plot
  		axis(1, at=1:4, labels=c("L3","6hrs", "24hrs","48hrs"), cex.axis=0.5, las=2)
  # 		goi <- c("maker-Draco-snap-gene-0.2143", "maker-Andromeda-augustus-gene-1.9563", "maker-Contig_17-augustus-gene-0.304")
  # 		for(j in 1:nrow(data_new)) {
  # 			lines(x=1:4, y=data_new[j,2:5], type="l", lwd=0.5, col="grey")
  # 		if(goi[1] %in% data_new$genes) #|| goi[2] %in% data_new$genes)
		# {	
		# 	lines(x=1:4, y=data_new[data_new$genes %in% goi[1],2:5], lwd=0.5, type="l", col="black")
		# 	#lines(x=1:4, y=data_new[data_new$genes %in% goi[2],2:5], lwd=0.5, type="l", col="brown")
			
		# }
		# else if(goi[2] %in% data_new$genes || goi[3] %in% data_new$genes)
		# {
		#  	lines(x=1:4, y=data_new[data_new$genes %in% goi[2],2:5], lwd=0.5, type="l", col="brown")
		#  	lines(x=1:4, y=data_new[data_new$genes %in% goi[3],2:5], lwd=0.5, type="l", col="green")	
		# }
		# else if(goi[3] %in% data_new$genes)
		# {
		# 	lines(x=1:4, y=data_new[data_new$genes %in% goi[3],2:5], lwd=0.5, type="l", col="green")
		# }
  		#}
  		par(new=TRUE)
  		plot(x=1:4, y=lambda_Djump[,i], pch=16, type = "b", axes = FALSE, bty = "n", xlab = "", ylab = "",
  			cex.lab=1, cex.axis=0.5, cex.main=1, col=colors[i])
		axis(side=4, at = zrange, las=2, cex.axis=0.5, cex.lab=1)
		mtext("lambda", side=4, line=3)
		
	}
	dev.off()
}