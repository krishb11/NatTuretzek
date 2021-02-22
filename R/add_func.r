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