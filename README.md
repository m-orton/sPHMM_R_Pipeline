## A kmer-Based Profile Hidden Markov Model (kPHMM) approach to Taxonomic Identification of DNA Barcode Data.

### To run this pipeline, the following files are needed in the same working directory in R/RStudio (R version 4.4.1 or higher):
kPHMM_Train&CV.R, kPHMM_Functions.R, seqstokmers.cpp and the folder containing sequence files in .fas or .fa format (zipped seqData folder provided for testing).

#### *** If running the pipeline to completion, ensure roughly 1 Gb of disk space is available per ~25k total sequences in your working directory.

### Package requirements for core functionality of the pipeline:
```
# BiocManager for Bioconductor packages if not installed
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install (if not installed) and load any bioconductor packages 
# (Biostrings is the only ones being used here)
# Define a function to check and install Bioconductor packages
load_bioPackages <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    BiocManager::install(package_name)
    library(package_name, character.only = TRUE)
  }
}

# Load the required Bioconductor packages
load_bioPackages("Biostrings")
load_bioPackages("treeio")

# List of CRAN packages used in this pipeline
packages <- c("DBI", "ape", "kmer", "aphid", "foreach", "readr", 
              "tidyverse", "plyr", "data.table", "splitTools", 
              "crayon", "future", "doFuture", "doRNG", "progressr", 
              "Rcpp", "cutpointr", "matrixTests", "rstatix", 
              "gtools", "RSQLite", "qs", "blaster","pryr", 
              "benchmarkme", "ggplot2", "gridExtra", "ggtree", 
              "tidytree", "taxonomizr", "colorspace","scales", 
              "magick", "viridis", "ggpubr")

# Install (if not installed) and load CRAN packages
for (i in 1:length(packages)) {
  if (!require(packages[i], character.only = TRUE)) {
    install.packages(packages[i], dependencies = TRUE)
    library(packages[i], character.only = TRUE)
  }
}

# One var for all packages used
packages <- sort(append(c("Biostrings", "treeio"), packages))

# Functions get loaded from a separate R script and must be loaded before the analysis proceeds
source("kPHMM_Functions.R")
```
