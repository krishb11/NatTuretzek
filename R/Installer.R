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
