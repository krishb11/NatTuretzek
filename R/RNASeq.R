library("sleuth")
library("biomaRt")



MyBat <- function(samples, pairwise_compare=TRUE, kal_dir, conditions, species) {
	value <- list(samples=samples, base_dir=kal_dir, conditions=conditions, species=species)
	if(pairwise_compare) {
		value$significant <- pairwiseCompare(value) 
	} else {
		value$significant <- GeneralCompare(value)
	}
	return(value)
}
pairwiseCompare <- function(BatObject) {
	mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = 'ensembl.org')
	t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart, useCache=FALSE)
	t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ext_gene = external_gene_name, ens_gene= ensembl_gene_id)
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

	#write results to csv
	sleuth_table <- sleuth_results(so, test = paste("conditions", BatObject$species, sep=""), show_all = TRUE)
	sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
	plot_pca(so, units = "est_counts", color_by = "conditions")
	plot_scatter(so)
	plot_vars(so)
	sleuth_live(so)
	return(sleuth_significant)
}

GeneralCompare <- function(BatObject) {
	return(NULL)
}


#data <- MyBat(samples=samples, kal_dir=kdir, conditions=c("mel", "mel","mel","suz", "suz", "suz", "bia", "bia", "bia"), species="suz")