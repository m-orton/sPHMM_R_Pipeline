# A kmer-Based Profile Hidden Markov Model (kPHMM) approach to Taxonomic Identification of DNA Barcode Data

###### Packages ######

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

# For performance benchmarking Tests
# install.packages('microbenchmark')
library(microbenchmark)

###### Warning/Message Settings ######

# To suppress summarize info in console for the dplyr package
options(dplyr.summarise.inform = FALSE)

# Max console output length for warnings - an expanded length is recommended for custom generated warnings in the console
options(warning.length = 8000L)

# To silence all warnings globally (not recommended since some of the functions below use warnings)
# options(warn=-1)

# Alternatively to turn warnings back on globally if disabled
# options(warn=0)

###### Working Directory ######

# ***Working directory must include these four files: PHMM_Train&Cross-Validate.R, PHHM_Test.R, PHMM_Functions.R, kmercntseq.cpp
# ***Working directory must also include a folder containing the fasta files to be run on this pipeline ex: seqData/ 
# ex: setwd("C:/PHMM_testrun)
setwd("~/")

##### Function Loading #####

# Functions get loaded from a separate R script and must be loaded before the analysis proceeds
source("PHMM_Functions.R")

##### Data Import #######

# First name a class to be used for training/testing from the provided fasta file, 
# ex: classRank <- "Arthropoda"

# ***This pipeline currently assumes only one class level taxonomy is being run at one time
classRank <- "Mammalia"

# ***The sequence data that is imported must be in .fa or .fas format
# ***The taxonomic name data labeling in the header info must be in this ordering:
# RecordID/SeqID must be first in the header info, then anywhere after in the header info Class;Order;Family;Genus;Species
# separated by either semi-colon, colon, underscore, space, tab or comma.
# example header info compatible with this pipeline: 
# RecordID;...;Class;Order;Family;Genus;Species OR RecordID;Class;Order;Family;Genus;Species OR RecordID;...;Class;Order;Family;Genus;Species;...

# fastaDir - the name of the folder where the fasta file(s) are located within the current working directory ex: seqData/
fastaDir <- "seqData/"

# Param 1 - fastaDir - name of the folder where the fasta file(s) are located
# Param 2 - classRank - character for the name of the class being run
dfTaxa <- fasta_Load(fastaDir, classRank)

###### Taxonomic Rank Selection ######

# If running through all four ranks set to All as the first param, 
# Alternatively if only one rank needs to be run, set subclassRank
# to ONE of the following: Order, Family, Genus or Species
# ex: subclassRank <- subclass_Rank("Species", classRank)

# Param 1 - ONE of the following: "All", "Order", "Family", "Genus" or "Species"
# Param 2 - character for the name of the class being run
subclassRank <- subclass_Rank("All", classRank)

###### Filtering by Sequence Quality ######

# Sequences can also be filtered for high gap content and/or high N content

# Threshold for proportion of N's present in each sequence - Default of 1% (0.01)
nFilter <- 0.01
# Threshold for proportion of gaps present in each sequence - Default of 1% (0.01)
gFilter <- 0.01

# Param 1 - dfTaxa - loaded sequence dataset 
# Param 2 - nFilter - numeric that represents the proportion of N's present in a sequence, if higher than this threshold, these sequences will be filtered out
# Param 3 - gFilter - numeric that represents the proportion of gap's present in a sequence, if higher than this threshold, these sequences will be filtered out
dfTaxa <- seq_Filter(dfTaxa, nFilter, gFilter)

###### Filtering by Taxonomic Rank ######

# Filter out groups with less sequences than the number specified by the groupMinNum param
# Default set to a minimum number of 25 sequences per taxonomic group, 
# ***This value MUST be set to a minimum of 2 before proceeding further since groups with only 1 sequence cannot 
# be used in downstream analysis of kmers (ttests to determine taxon-specificity of kmers require more than 1 observation/sequence)
groupMinNum <- 25

# In addition, to limit the size of datasets and speed up downstream kmer analyses, the
# groupMaxNum param is used - Default set to a maximum number of 50 sequences per taxonomic group, 
# ***Sequences are randomly sampled (n=groupMaxNum) for each taxonomic group IF a group has a greater number of sequences than the groupMaxNum value set 
# ***groupMaxNum can also be set to NA if an upper limit to the number of sequences is not needed (ex: relatively small datasets of less than 1000 sequences)
groupMaxNum <- 50

# Set a global seed for any commands that use random sampling 
set.seed(1)

# Param 1 - dfTaxa - loaded sequence dataset 
# Param 2 - subclassRank - taxonomic ranks being analyzed 
# Param 3 - groupMinNum - minimum number of sequences taxonomic groups must have to be included in the analyses
# Param 4 - groupMaxNum - maximum number of sequences taxonomic groups can have before being randomly sampled (n=the value specified in this param)
dfTaxa <- subclass_Filter(dfTaxa, subclassRank, groupMinNum, groupMaxNum)

###### Parallel Processing ######

# Set to TRUE to enable parallel processing for wilcoxon rank-sum tests (discriminitive kmer determination) and PHHM setup/training 
# ***Recommended for faster runtimes but comes at the cost of higher overall memory usage
runParallel <- TRUE

# numCores param - Default set to available cores - 1 if runParallel == TRUE
# If set to FALSE, numCores will always be set to 1
# ***Recommended to always subtract 1 core if runParallel == TRUE and this program is being run locally (so that 1 CPU is available for local applications)
if(runParallel == TRUE){
  numCores <- availableCores()
} else {
  numCores <- 1
}

###### Cross-Validation Setup ######

# Creating folds for cross-validation. Subdividing into folds is intended for cross-validation tests
# ***If numFolds is set to 1, folds will not be created and cross-validation procedures will not occur

# numFolds - param for total number of folds to create, default set to 5 folds, set to 1 for no fold creation
# and numFolds MUST be set to a minimum of 3 if performing cross-validation tests, thus numFolds cannot equal 2
numFolds <- 5

# testNum, trainNum - params for what proportion of folds to be used for training and testing
# These params must be integers summing to numFolds. These params are only relevant for cross-validation and ignored otherwise
# ex: if numFolds = 5: default proportions of test/train are testNum = 2, trainNum = 3
testNum <- 2
trainNum <- 3

# Set seed for each fold in the stratified sampling process when performing the cross-validation step in the stratify_Sampling function
# This is ignored if numFolds = 1
crossValSeed <- 1:numFolds

# Stratified Sampling of dfTaxa dataframe - if not performing cross-validation tests, dfTaxa will not be assigned folds

# For stratified sampling, the following params are needed:
# Param 1 - dfTaxa - sequence dataset being used
# Param 2 - classRank - character for the class-level taxonomy currently being run
# Param 3 - subclassRank - character for the all subclass ranks being run
# Param 4 - crossValSeed - numeric vector for the seed set to each fold in cross-validation
# Params 5-7 - numeric vectors for cross-validation params - numFolds, testNum & trainNum (only used if numFolds >= 3)
stratify_Sampling(dfTaxa, classRank, subclassRank, crossValSeed, numFolds, testNum, trainNum)

###### Kmer Extraction ######

# Kmer length param to specify which kmer lengths to run
klen <- c(2,3,4)

# Iteration over each kmer length with the kmer_Find function for kmer lengths greater than 1

# For kmer extraction, the following params are needed:
# Param 1 - dfTaxa - sequence dataset being used
# Param 2 - dir - character for the main directory where the analysis is being performed
# Param 3 - klen - numeric vector of specified kmer length(s) greater than 1
kmerFiles <- kmer_Find(dfTaxa, dir, klen)

# *** The kmers extracted using the kmer_Find function are written to the disk in their own directory depending on taxonomic rank, 
# they will be located in the KmerData subdirectory. 
# ***The final output of this function is a list consisting of both the raw kmer filepaths and the chunked kmer filepaths

###### Discriminitive Kmer Determination ######

# For kmer lengths greater than 1 are being run, run the kmer_discrimFind function with the params specified below

# Param for filtering by p-value significance for the wilcoxon rank sum test (Default = 0.01).
sigVal <- 0.01

# Param for the minimum number of discriminative kmers that must exist for each group and taxonomic rank, 
# min value for the minimum number must be 1 (Default = 1).
numMinK <- 1
  
# Param for the maximum number of discriminitive kmers allowed for each group and taxonomic rank, 
# discrim kmers will first be sorted by significance and numMaxK will govern how many of the kmers (with the highest significance) to choose from once sorted
# Can alternatively set to NA in which case this param will be ignored (Default = 25).
numMaxK <- 25
  
# Params for kmer_discrimkFind:
# Param 1 - kmerFiles - list of filenames corresponding to the directory where the .tsv files for each raw kmer dataset are located
# Param 2 - dfTaxa - dataframe containing all sequence data
# Param 3 - numFolds - numeric for specified number of folds used
# Param 4 - subclassRank - character for taxonomic ranks to be run
# Param 5 - klen - numeric for specified kmer length(s)
# Param 6 - sigVal - numeric for the sig threshold set for wilcox rank sumtests
# Param 7 - numMinK - numeric for the min number of kmers allowed for each taxonomic group at each rank and klen
# Param 8 - numMaxK - numeric for the max number of kmers allowed for each taxonomic group at each rank and klen
# Param 9 - runParallel - boolean indicating whether to use parallel processing for ttests 
# Param 10 - numCores - numeric indicating the number of cores to use
  
# Finally we run through each kmer length and taxonomic rank to find discriminitive kmers for each group, 
# the kmer_discrimkFind function will load in chunks one by one until a minimum number (numMinK) of discriminitive kmers can be found for each group. 
kmerResults <- kmer_discrimFind(kmerFiles, dfTaxa, numFolds, subclassRank, klen, 
                                sigVal, numMinK, numMaxK, runParallel, numCores)
  
# ***The final output of this function is a list consisting of two elements: 
# Results - dataframe detailing the discriminitive kmers discovered from the statistical testing performed with their associated p-values
# Tally - dataframe consisting of tallied counts of discriminitive kmers found for each group and rank

###### kPHMM Setup & Training ###### 

# numMaxK param - maximum number of discrim kmer positions found per group/rank/klen that are used in 
# kPHHM model setup (after positions are ordered by highest discrim kmer count), Default = 200
# This is mainly to limit the complexity of each kPHMM model (as some models may have several hundreds of potential discrim kmer positions to work with) 
# while also reducing the runtime, memory and diskspace usage during the setup/training process
numMaxKPos <- 200
  
# numMinKSeq param - minimum number of discrim kmers that must be present in each sequence before model setup, Default = 1
# This param is necessary due to the fact that not every sequence belonging to a certain group will necessarily possess the discrim kmers identified for that particular group
numMinPropKSeq <- 1
  
# gapThreshold param - percentage (as a decimal) of acceptable gap content across all sequences at each discrim kmer position during model setup, Default = 0.8
# Higher values result in increased tolerance for gaps (ex: 0.9 means 90% or lower gap content per kmer position is acceptable), 
# Lower values result in decreased tolerance for gaps (ex: 0.2 means 20% or lower gap content per kmer position is acceptable)
gapThresholdKPos <- 0.8

# numMaxiter param - maximum number of EM iterations the model must perform before the cycling process is terminated and the model is returned. Default = 100
numMaxiter <- 100

# deltaLL param - represents as a numeric the maximum change in log likelihood between EM iterations (during kPHMM model training) before
# model convergence is attained. Defaults to 1E-07. 
deltaLL <- 1E-07
  
# Params for kphmm_Train:
# Param 1 - kmerResults - list of discrim kmer results
# Param 2 - kmerFiles - list of filenames corresponding to the directory where the .tsv files for each raw kmer dataset are located
# Param 3 - klen - numeric for specified kmer length(s)
# Param 4 - numFolds - numeric for specified number of folds used
# Param 5 - chunkSize - numeric for the total number of rows from the raw kmer dataset to be read into memory in one chunk
# Param 6 - numMaxKPos - numeric for the max number of discrim kmer positions that are used during model setup
# Param 7 - mumMinKSeq - numeric for the min number of discrim kmers that must be present in each sequence during model setup
# Param 8 - gapThreshold - numeric for the percentage (as a decimal) of acceptable gap content across all sequences at each discrim kmer position during model setup
# Param 9 - numMaxiter - numeric for the maximum number of EM iterations the model must perform before the cycling process is terminated and the model is returned 
# Param 10 - subclassRank - character for taxonomic ranks to be run
# Param 11 - deltaLL - numeric for the the maximum change in log likelihood between EM iterations before the cycling procedure is terminated
# Param 12 - numCores - numeric indicating the number of cores to use if runParallel is TRUE, if FALSE, numcores will always equal 1
  
# Model setup and training of PHHM's using all of the specified parameters above
modelList <- kphmm_Train(kmerResults, kmerFiles, klen, numFolds, 
                         chunkSize, numMaxKPos, numMinKSeq, gapThreshold,
                         numMaxiter, subclassRank, deltaLL, numCores)
  
# ***Upon completion, trained models will then get written to the ModelData/ folder in .hmm format
# ***The final output of this function is a list consisting of three elements: 
# Model Overview - dataframe indicating which groups successfully had kPHMM models trained for them 
# Model Data - Model kphmms for each group including the filepaths of each model and input matrix written to disk
# Model Input Matrices - The input matrices used for each group during model setup

###### kPHMM Testing ######

# ***The ClassifyCV function will be used when performing cross-validation tests, otherwise the kphmm_ClassifyI function should be used when not performing cross-validation tests,

# if klengths > 1 are being run, this function will classify all "test" labeled sequences for each fold
# Params for kphmm_ClassifyCV:
# Param 1 - kmerFiles - list of filenames corresponding to the directory where the .tsv files for each raw kmer dataset are located
# Param 2 - modelList - list of dataframes and filenames corresponding to each of the kphmm models trained
# Param 3 - numFolds - numeric for specified number of folds used
# Param 4 - subclassRank - character for taxonomic ranks to be run
# Param 5 - klen - numeric for specified kmer length(s)
# Param 6 - classRank - character for class level taxonomy being run
  
# Test each fold combination for each rank on the kPHMM models that have been set up and trained
modelResults <- kphmm_ClassifyCV(kmerFiles, modelList, numFolds, subclassRank, klen, classRank)
  
# ***The final output of this function is a list consisting of 2 elements: 
# Model Results Overview - dataframe indicating the % correct ids per group for each subclassRank
# Model Plot filepaths - list consisting of NULL values until plotResults function is run in which case the 
# filepaths of results plot written to disk will be included in the modelResults

###### Plotting Results - Optional ######

# Data Visualization
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

# Function for plotting
source("PHMM_plotting_functions.R")

# plotResults function - if run, will plot a heatmap (% correct ids), ROC curve, MCC curve and F1 curve (across folds) of the results data for each taxonomic rank and klen, 
# write each one to a .png file and store them all in a subfolder labeled ModelResults

# Numplots is an estimate of how many plots to subdivide the %correct ids results by on a per rank basis depending on the size of the dataset
# ex: numPlots <- c(1, 2, 6, 16) would mean 1 %correct id plot for Order, 2 for Family, 6 for Genus, 16 for Species
numPlots <- c(1, 2, 6, 16)

# Param 1 - modelResults list generated from kphmm_ClassifyCV function
# Param 2 - numPlots numeric for a rough estimate of how many plots per rank
# Param 3 - klen numeric for the specified kmer lengths

# Run function for plotting results
modelResults <- plot_Results(modelResults, numPlots, klen)

# ***The final output of this function is a list consisting of 2 elements: 
# Model Results Overview - dataframe indicating the % correct ids per group for each subclassRank (unchanged from previous function)
# Model Plot filepaths - the filepaths of results plots written to disk (in .png format) will be included in the modelResults
