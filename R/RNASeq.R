library("sleuth")
library("biomaRt")
library("dplyr")


#' Constructor and automated function to perform pairwise comparison of RNASeq data.
#' æparam samples: A character vector with all the kallisto directories of interest in order
#' @param kal_dir: path to the base directory where all the kallisto directories are present
#' @param conditions: A vector that defines the condition of each sample
#' @param species: A vector defining he species from which the samples are derived in a one-to-one manner
#' @retuns: A list with all the input parameters and the output of downstream RNAseq analysis in sleuth that can be viewed using sleuth_live function.
MyBat <- function(samples, pairwise_compare=TRUE, kal_dir, conditions, species) {
	value <- list(samples=samples, base_dir=kal_dir, conditions=conditions, species=species, significant=NULL)
	if(pairwise_compare) {
		value$significant <- pairwiseCompare(value) 
	} else {
		value$significant <- GeneralCompare(value)
	}
	return(value)
}
pairwiseCompare <- function(BatObject) {
	t2g <- read.csv("/home/krishna/test_Nat/run_Kallisto/t2g.csv")
	t2g <- dplyr::rename(t2g, target_id = genes, ext_gene = names)
	#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = 'ensembl.org')
	#t2g <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)
	#t2g <- dplyr::rename(t2g, target_id = ensembl_gene_id, ext_gene = external_gene_name)
	s2c <- data.frame(path=as.character(BatObject$base_dir), sample=as.character(BatObject$samples), conditions=as.character(BatObject$conditions))
	s2c$path <- as.character(s2c$path)
	s2c$sample <- as.character(s2c$sample)
	s2c$conditions <- as.character(s2c$conditions)
	print(s2c)
	#build sleuth object
	so <- sleuth_prep(s2c, target_mapping = t2g, extra_bootstrap_summary = TRUE, transformation_function = function(x) log2(x + 0.5), full_model = ~conditions, read_bootstrap_tpm = TRUE)
	#test for significance
	so <- sleuth_fit(so)
	models(so)
	so <- sleuth_wt(so, which_beta = paste("conditions", BatObject$species, sep=""))
	so <- sleuth_fit(so, ~1, "reduced")
	so <- sleuth_lrt(so, 'reduced', 'full')
	results_wt <- sleuth_results(so, test = paste("conditions", BatObject$species, sep=""), show_all = TRUE)
	results_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
	lrt.sig_ids <- results_lrt$target_id[which(results_lrt$qval < 0.05)]
	wt.sig_ids <- results_wt$target_id[which(results_wt$qval < 0.05)]
	ids <- wt.sig_ids[wt.sig_ids %in% lrt.sig_ids]
	head(ids)
	shared_results <- results_wt[results_wt$target_id %in% ids,]
	write.csv(shared_results, file=BatObject$output)
	write.csv(results_wt, file=paste("WT_", BatObject$output, sep=""))

	pdf(BatObject$output)
	plot_pca(so, units = "est_counts", color_by = "conditions")
	plot_scatter(so)
	plot_vars(so)
	dev.off()
	return(so)
}

GeneralCompare <- function(BatObject) {
	return(NULL)
}


#data <- MyBat(samples=samples, kal_dir=k, conditions=c("mel", "mel","mel", "bia", "bia","bia"), species="mel", test="bia_48_v_mel.csv")
#data <- MyBat(samples=samples, kal_dir=k, conditions=c("female", "female", "female", "male", "male","male"), species="male", test="48_v_female.csv")	
#data <- MyBat(samples=samples, kal_dir=k, conditions=c("suz", "suz", "suz", ), species="suz", test="48_v_mel.csv")

#samples <- c( "bia_24_female1", "bia_24_female2", "bia_24_female3", "bia_48_female1" ,"bia_48_female2", "bia_48_female3","bia_6_female1" ,"bia_6_female2", "bia_6_female3", "bia_L3_female1", "bia_L3_female2", "bia_L3_female3",
#	"mel_24_female1", "mel_24_female2", "mel_24_female3","mel_48_female1" ,"mel_48_female2", "mel_48_female3", "mel_6_female2", "mel_6_female3", "mel_L3_female1", "mel_L3_female2", "mel_L3_female3", 
# "suz_24_female1", "suz_24_female2", "suz_24_female3", "suz_48_female1" ,"suz_48_female2", "suz_48_female3","suz_6_female1" ,"suz_6_female2", "suz_6_female3", "suz_L3_female1", "suz_L3_female2", "suz_L3_female3")

#species <- c(rep(c("bia"), each=12), rep(c("mel"), each=11), rep(c("suz"), each=12))
#conditions("48hrsAPF", "24hrsAPF", "6hrs")



#mel vs bia - check with suz
#up and downregulated genes and compare with deseq and sleuth
#Visual memory of gential disc. repression should be explained clearly
#ac,sc,wg - highest in mel, med in bia, lowest in suz
#dpp high in biarmipes only 
#Venn diagrams in all comparsions with up and downregulation separately.
#networks if possible

#conditions <- c( rep(c("24hrsAPF"), each=3), rep(c("48hrsAPF"), each=3), rep(c("6hrsAPF"), each=3), rep(c("L3"), each=3),
#rep(c("24hrsAPF"), each=3), rep(c("48hrsAPF"), each=3), rep(c("6hrsAPF"), each=2), rep(c("L3"), each=3),
#rep(c("24hrsAPF"), each=3), rep(c("48hrsAPF"), each=3), rep(c("6hrsAPF"), each=3), rep(c("L3"), each=3))




#read text files and csv files

