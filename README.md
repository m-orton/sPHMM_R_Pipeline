A kmer-Based Profile Hidden Markov Model (kPHMM) approach to Taxonomic Identification of DNA Barcode Data.

To run this pipeline, the following files are needed in the same working directory in R/RStudio (R version 4.2.1 or higher):
kPHMM_Train&Test.R 
PHMM_Functions.R 
kmercountwritecolumns.cpp 
A folder containing sequence files in .fas or .fa format (zipped seqData folder provided for testing).

Optionally, for plotting of tests results from the pipeline, the PHMM_plotting_function.R file is required.

Package requirements for core functionality of the pipeline:

# Bioinformatics and Sequence Handling
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)
# BiocManager::install("GenomicRanges", force = TRUE)
library(GenomicRanges)
# DNAbin format
# install.packages("ape")
library(ape)

# Functional Programming and Iteration
# install.packages("foreach")
library(foreach)

# File writing in tsv/csv and 
# reading in chunked data files
# install.packages("readr")
library(readr)
# install.packages("chunked)
library(chunked)

# Efficient Data Manipulation
# install.packages("tidyverse")
library(tidyverse)
# install.packages("plyr")
library(plyr)
# install.packages("dplyr")
library(dplyr)
# install.packages("data.table")
library(data.table)
# install.packages("tibble")
library(tibble)
# install.packages("purr")
library(purrr)
# install.packages("stringr")
library(stringr)
# install.packages("tidyr")
library(tidyr)

# Cross-validation setup
# install.packages("splitTools")
library(splitTools)

# Parallel processing functions
# install.packages("future")
library(future)
# install.packages("doFuture")
library(doFuture)
# install.packages('doRNG')
library(doRNG)

# Kmer count matrices
# install.packages("kmer")
library(kmer)
# Faster matrix col sums
# install.packages("Rfast")
library(Rfast)
# C++ Integration
# install.packages("Rcpp")
library("Rcpp")

# For performing fast statistical testing on matrices of kmer data
# install.packages("matrixTests")
library(matrixTests)

# Pvalue adjustment
# install.packages("rstatix")
library(rstatix)

# Sequence Alignment and Profile HMMs
# install.packages("aphid")
library(aphid)

# Kmer Permutations
# install.packages("gtools")
# library(RcppAlgos)
library(gtools)

# Custom coloring of messages with cat
# install.packages("crayon")
library(crayon)

For plotting, additional packages are required:

# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggforce")
library(ggforce)
library(lattice)
# install.packages("grid")
library(grid)
# install.packages("gridExtra")
library(gridExtra)
# install.packages("gridtext")
library(gridtext)
# install.packages("viridis")
library(viridis)
# install.packages("cutpointr")
library(cutpointr)
# install.packages("gtable")
library(gtable)
# install.packages("egg")
library(egg)
# install.packages("cowplot")
library(cowplot)
