#Install all additional necessary packages

packages = c("dplyr", "scales", "biomaRt",
             "sleuth", "HTSFilter", "HTSCluster", "Biostrings",
              "DESeq2", "rlist", "readr", "ggplot2")

## Now load or install&load all
packages.list <- lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    } else {
    	library(x, character.only=TRUE)
    }
  })

system("conda install kallisto")



# Monday 10'o clock thesis
# check the forceps and trash useless 
# Leftovers leave on a separate rack
# clean up th boxes in teh freezer and fridge
# main clones naming and excel sheet updates

