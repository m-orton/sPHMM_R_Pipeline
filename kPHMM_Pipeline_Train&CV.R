# A kmer-Based Profile Hidden Markov Model (kPHMM) approach to taxonomic identification of DNA sequence data - Training & Cross-validation 

# *** If running the pipeline to completion, ensure roughly 1 Gb of disk space is available per ~25k total sequences in your working directory. ***

###### 1 - Packages & Functions for Pipeline Functionality ######

# *** Ensure R version is updated to 4.4.1 prior to using this pipeline ***

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

###### 2 - Warning/Message Settings ######

# To suppress summarize info in console for the dplyr package
options(dplyr.summarise.inform = FALSE)

# To silence all warnings globally (not recommended since some of the functions below use warnings)
# options(warn=-1)

# Alternatively to turn warnings back on globally if disabled
# options(warn=0)

###### 3 - Global Seed & Important Params ######

# Set a global seed for any general commands that use random sampling 
set.seed(1)

# *** The following params below are being used primarily for testing purposes and can be disabled otherwise. ***

# klen1Compare - if TRUE, will also train standard PHMM models on each group and rank to make comparisons in accuracy and performance
# when compared to the kmer-based models
klen1Compare <- TRUE

# runBenchmark - if TRUE, will record memory and processing times for each major function, note that performing benchmarking will increase 
# pipeline runtime to a minor degree
runBenchmark <- TRUE

###### 4 - Working Directory ######

# *** Working directory must include these three files: kPHMM_Pipeline_Train&CV.R, kPHMM_Functions.R, seqstokmers.cpp
# and must also include a folder containing the fasta files to be run on this pipeline ex: seqData/
# ex: setwd("C:/kPHMM_testrun) ***
setwd("D:/Chap3&4/Aug5Scripts")

##### 5 - Data Import #######

# First name a class to be used for training/testing from the provided fasta file, 
# ex: classRank <- "Arthropoda"

# *** This pipeline currently assumes only one class level taxonomy is being run at one time ***
classRank <- "Amphibia"

# *** The sequence data that is imported must be in .fa or .fas format
# and the taxonomic name data labeling in the header info must be in this ordering:
# RecordID/SeqID must be first in the header info, then anywhere after in the header info Class;Order;Family;Genus;Species
# separated by either semi-colon, colon, underscore, space, tab or comma.
# examples of header info compatible with this pipeline with different separator chars: 
# RecordID;...;Class;Order;Family;Genus;Species OR RecordID;Class;Order;Family;Genus;Species OR RecordID;...;Class;Order;Family;Genus;Species;... 
# RecordID,...,Class,Order,Family,Genus,Species OR RecordID,Class,Order,Family,Genus,Species OR RecordID,...,Class,Order,Family,Genus,Species,... ***

# fastaDir - the name of the folder where the fasta file(s) are located within the current working directory ex: seqData/
fastaDir <- "seqData/"

# Param 1 - fastaDir - name of the folder where the fasta file(s) are located
# Param 2 - classRank - character for the name of the class being run
dfTaxa <- fasta_Load(fastaDir, classRank)

###### 6 - Taxonomic Rank Selection ######

# If running through all four ranks set to All as the first param, 
# Alternatively if only one rank needs to be run, set subclassRank
# to ONE of the following: Order, Family, Genus or Species
# ex: subclassRank <- subclass_Rank("Species", classRank)

# Param 1 - ONE of the following: "All", "Order", "Family", "Genus" or "Species"
# Param 2 - character for the name of the class being run
subclassRank <- subclass_Rank("All", classRank)

# *** Note that running the pipeline at one rank only will limit some options in the downstream analysis:
# Classification strategy options are limited only to the MaxScore option and this may impact classification accuracy
# Results plotting options are more limited - see section 15 - Results Plotting

###### 7 - Filtering by Sequence Length and Quality ######

# *** Currently only testing pipeline with approx. standard length barcode sequences of 600-700 bp
# but will be expanded to variable length sequences of greater difference in length soon and does technically 
# work with variable length sequences. ***

# Threshold for proportion of N's present in each sequence - Default of 1% (0.01) of total sequence length or less. 
# Set to NA to ignore
nFilter <- 0.01
# Threshold for proportion of gaps present in each sequence - Default of 1% (0.01) of total sequence length or less.
# Set to NA to ignore
gFilter <- 0.01

# Sequence length range to consider for train datasets (and test datasets if performing cross-validation), 
# Set both values to NA to ignore
seqLenRange <- c(600, 700)

# Param 1 - dfTaxa - loaded sequence dataset 
# Param 2 - nFilter - numeric that represents the proportion of N's present in a sequence, if higher than this threshold, these sequences will be filtered out
# Param 3 - gFilter - numeric that represents the proportion of gap's present in a sequence, if higher than this threshold, these sequences will be filtered out
# Param 4 - seqLenRange -  numeric for the sequence length range to consider for train datasets (and test datasets if performing cross-validation)
dfTaxa <- seq_Filter(dfTaxa, nFilter, gFilter, seqLenRange)

###### 8 - Filtering by Number of Sequences per Taxonomic Rank ######

# Filter out groups with less sequences than the number specified by the groupMinNum param
# Default set to a minimum number of 10 sequences per taxonomic group, 
# *** This value MUST be set to a minimum of 2 before proceeding further since groups with only 1 sequence cannot
# be used in downstream analysis of kmers (wilcoxon rank sum tests to determine taxon-specificity of kmers require 
# more than 1 observation/sequence). ***
groupMinNum <- 5

# In addition, to limit the size of datasets and speed up downstream kmer analyses, the groupMaxNum param is used - Default set to a maximum number of 100 sequences 
# per taxonomic group. This value applies to species level if running all taxonomic ranks, thus groupMaxNum for order for example would be groupMaxNum * number of orders * 
# number of families * number of genuses whereas groupMaxNum for family would be groupMaxNum * number of families * number of genuses and so on.
# *** Sequences are randomly sampled (n=groupMaxNum) for each taxonomic group IF a group has a greater number of sequences than the groupMaxNum value set.
# groupMaxNum can also be set to NA if an upper limit to the number of sequences is not needed. ***
groupMaxNum <- NA

# Param 1 - dfTaxa - loaded sequence dataset 
# Param 2 - subclassRank - taxonomic ranks being analyzed 
# Param 3 - groupMinNum - the minimum number of sequences that a taxonomic group must have in order to be included in the analysis
# Param 4 - groupMaxNum - the maximum number of sequences that a taxonomic group can have before being randomly sampled (n=the value specified in this param)
dfTaxa <- subclass_Filter(dfTaxa, subclassRank, groupMinNum, groupMaxNum)

###### 9 - Cross-Validation Setup ######

# Creating folds for cross-validation. Subdividing into folds is intended for cross-validation tests
# *** If numFolds is set to 1, folds will not be created and cross-validation procedures will not occur downstream in the analysis ***

# numFolds - param for total number of folds to create, default set to 5 folds, set to 1 for no fold creation
# and numFolds must be set to a minimum of 3 if performing cross-validation tests
numFolds <- 5

# testNum, trainNum - params for what proportion of folds to be used for training and testing
# These params must be integers summing to numFolds. These params are only relevant for cross-validation and ignored otherwise
# ex: if numFolds = 5: default proportions of test/train are testNum = 2, trainNum = 3
testNum <- 2
trainNum <- 3

# Set a seed for each fold in the stratified sampling process when performing the cross-validation step in the stratify_Sampling function
# This param is ignored if numFolds = 1
crossValSeed <- 1:numFolds

# Stratified Sampling of dfTaxa dataframe - if not performing cross-validation tests, dfTaxa will not be assigned folds

# For stratified sampling, the following params are needed:
# Param 1 - dfTaxa - sequence dataset being used
# Param 2 - classRank - character for the class-level taxonomy currently being run
# Param 3 - subclassRank - character for the all subclass ranks being run
# Param 4 - crossValSeed - numeric vector for the seed set to each fold in cross-validation
# Params 5-7 - numeric vectors for cross-validation params - numFolds, testNum & trainNum (only used if numFolds >= 3)
stratify_Sampling(dfTaxa, classRank, subclassRank, crossValSeed, 
                  numFolds, testNum, trainNum)

###### 10 - Parallel Processing ######

# Set to TRUE to enable parallel processing for discrim kmer finding/PHHM setup/training and testing/cross-validation.
# *** Recommended for faster runtimes but comes at the cost of higher overall memory usage. ***
runParallel <- TRUE

# if TRUE, set global settings for parallel processing
if(runParallel == TRUE){
  # Option for future R package that indicates the maximum RAM that can be allocated per core, 1.5 Gb seems to be sufficient in most cases
  # but may need to be set higher if the sequence dataset is quite large (>25k sequences)
  futureGlobals <- 1500*1024^2
  # Seed for parallel processing with doRNG 
  doRNGseed <- 1
  # numCores param - Default set to available cores if runParallel == TRUE
  # If set to FALSE, numCores will always be set to 1. Recommended to subtract 1 core if 
  # runParallel == TRUE and this program is being run locally
  numCores <- availableCores()
  # False, numCores is always set to 1
} else {
  numCores <- 1
  futureGlobals <- NA
  doRNGseed <- 1
}

###### 11 - Kmer Extraction ######

# Kmer length param to specify which kmer lengths to run,  maximum kmer length is 11
# It is recommended to provide a range of kmer lengths as this will likely improve model accuracy as some groups perform better with lower 
# kmer lengths while others perform better at higher kmer lengths, default kmer length range is 2-5.
klen <- c(2,3,4,5)

# Add klen 1 if klen1Compare is TRUE. In the case of klen 1, the kmer_Find function will simply convert the sequence dataset
# into a matrix-like data structure with each column representing the positions of nucleotides across the full length of the sequence - 
# this data structure is what is used in downstream analysis for klen 1.
if(klen1Compare == TRUE){
  klen <- append(1, klen)
}

# For kmer extraction, the following params are needed:
# Param 1 - dfTaxa - sequence dataset being used
# Param 2 - subclassRank - character for taxonomic ranks to be run
# Param 3 - dir - character for the name of db where kmer data is being written to
# Param 4 - klen - numeric vector of specified kmer length(s)
# Param 5 - runBenchmark - boolean indicating whether to run benchmarking tests (benchmark results also get written to db)
# Param 6 - runParallel - boolean indicating whether parallel processing should be used
# Param 7 - numCores - numeric indicating the number of cores to use if runParallel is TRUE, if FALSE, numcores will always equal 1
# Param 8 - futureGlobals - numeric indicating the maximum amount of RAM that is allocated per core when runParallel is TRUE
# Param 9 - doRNGseed - numeric for the seed set for dorng foreach command
kmerData <- kmer_Find(dfTaxa, subclassRank, dir, klen, runBenchmark, 
                      runParallel, numCores, futureGlobals, doRNGseed)

# The kmers with their associated frequencies and positions are written to a sqlLite db. 
# The final output of this function is a list consisting of both the kmer frequency
# db table names and the kmer position db table names.

# The kmer data from the db can be accessed with the following examples below:

# Examples below when not using cross-validation: 
# table <- read_dbTable(dbPath, "raw_klen_2_cnt") or table <- read_dbTable(dbPath, "raw_klen_3_pos")

# Examples when using cross-validation:
# table <- read_dbTable(dbPath, "raw_klen_2_cnt") or 
# table <- read_dbTable(dbPath, "raw_klen_3_f1_train_pos") or
# table <- read_dbTable(dbPath, "raw_klen_3_f1_test_pos")

# The db will also contain a table with benchmarking results if runBenchmark is TRUE.
# The benchmarking data from the db can be accessed easily: table <- read_dbTable(dbPath, "bench_kmer_find")

###### 12 - Discriminative/Taxon-Specific Kmer Determination ######

# Different klen var for klen without 1 since klen 1 is not used in this section
klenSub <- klen[klen!=1]

# Params for filtering by p-value significance for the wilcoxon rank sum test, Default = 0.05 for maxSigVal.
# The kmer_discrimFind function will find all discrim kmers at equal to or below the specified threshold (maxSigVal) and 
# the minimum number of discrim kmers that must be present for each group/fold/klen (otherwise they are filtered out) is set 
# by the numMinK param, Default = 1, min = 1. numMinK cannot be set higher than numMaxK for obvious reasons (see below).
maxSigVal <- 0.05
numMinK <- 1

# Param for the maximum number of discriminative kmers permitted for each group and taxonomic rank, 
# discrim kmers will be sorted by significance and kCoverage and numMaxK will govern how many of the top-ranked kmers 
# (with the highest pvalue significance and lowest kCoverage) to choose from once sorted. This value scales with
# klen, it is purposely set lower for lower klens and increases proportionally with increasing klen. This is because the number of possible
# kmer permutations increases exponentially with increasing klen. A scaling coefficient is used to determine numMaxK for each klen, 
# and is currently set to roughly 20% of all kmer permutations for klen=2 which means numMaxK=4, 15% for klen=3 which means numMaxK=10, 
# 10% for klen=4 which means numMaxK=26, 5% of kmers for klens=5 which means numMaxK=52, and 1% of kmers for klen>5. 
scaleCoefMaxK <- c(0.2, 0.15, 0.1, 0.05, 0.01)
numMaxK <- calc_numMaxK(klenSub, scaleCoefMaxK)

# Params for kmer_discrimkFind:
# Param 1 - kmerData - list of db table names corresponding to each kmer length for frequency and position data
# Param 2 - dfTaxa - dataframe containing all sequence data
# Param 3 - numFolds - numeric for specified number of folds used
# Param 4 - subclassRank - character for taxonomic ranks to be run
# Param 5 - klenSub - numeric for specified kmer length(s) without klen=1, min of klen=2, max of klen=11
# Param 6 - maxSigVal - numeric for the upper limit of acceptable pvalues for discrim kmers
# Param 7 - numMinK - numeric for the min number of kmers permitted for each taxonomic group at each rank and klen
# Param 8 - numMaxK - numeric for the max number of kmers permitted for each taxonomic group at each rank and klen
# Param 9 - runParallel - boolean indicating whether parallel processing should be used
# Param 10 - numCores - numeric indicating the number of cores to use if runParallel is TRUE, if FALSE, numcores will always equal 1
# Param 11 - futureGlobals - numeric indicating the maximum amount of RAM that is allocated per core when runParallel is TRUE
# Param 12 - dbPath - character for path to SQLite db
# Param 13 - runBenchmark - boolean indicating whether to run benchmarking tests

# Finally we run through each kmer length and taxonomic rank to find discriminitive kmers for each group, 
# the kmer_discrimkFind function will load in chunks one by one until a minimum number (numMinK) of discriminitive kmers can be found for each group. 
kmerResults <- kmer_discrimFind(kmerData, dfTaxa, numFolds, subclassRank, 
                                klenSub, maxSigVal, numMinK, numMaxK, runParallel, 
                                numCores, futureGlobals, dbPath, runBenchmark)

# The final output of this function is a list consisting of two elements: 
# Results - dataframe detailing the discriminative kmers discovered from the statistical testing performed with their associated p-values
# Tally - dataframe consisting of tallied counts of discriminative kmers found for each group and rank.
# Fail - dataframe consisting of groups at each rank, klen and fold (if applicable) where kmers could not be found
# BenchData written to db if runBenchmark is TRUE, this data can be accessed easily: table <- read_dbTable(dbPath, "bench_kmer_discrimFind")

# More info on the discrimFind function and kCoverage:
# kCoverage (as seen in the kmer results data) is a measure of kmer taxon-specificity that represents the proportion of groups at a certain taxonomic rank
# that possess a specific kmer after the wilcoxon test pvalue cutoffs are applied, ex: if a kmer is found to possess a kCoverage of 0.1 after its pvalue cutoff is applied, 
# this will mean that only 10% of all groups at a specific rank (and fold if performing cross-val tests) possess that kmer and the kmer_discrimFind function 
# is optimized for finding kmers that possess the highest significance in terms of pvalue and also possess the lowest possible kCoverage to maximize taxon-specificity.

###### 13 - kPHMM Setup & Training ######

# gapThreshold param - the acceptable proportion of gaps at each model position during training. Only positions that meet this threshold will 
# be counted as match states. Defaults to 0.5, range 0-1. Higher values mean less stringent tolerance for gaps at every position, 
# lower values mean more stringent tolerance for gaps at every position. A value of 0 would mean all positions that contain gaps 
# would not count as match states, a value of 1 would mean that all positions with gaps would be counted as match states. 
gapThreshold <- 0.5

# numMaxiter param - maximum number of EM iterations the model can perform before the cycling process is terminated and the model is returned. Default = 100, min = 1, 
# setting this value very high will increase time to convergence significantly (with eventual diminishing returns in model precision) but setting it too low will 
# likely yield less precision in model training.
numMaxiter <- 100

# deltaLL param for model setup - represents as a numeric the maximum change in log likelihood between EM iterations (during kPHMM model setup) before
# model convergence is attained. Defaults to 1E-07. This value is used to monitor the convergence of the training process. During training, 
# the goal is to maximize the log-likelihood. As training progresses, the changes in log-likelihood should become smaller and smaller indicating 
# that the model parameters are stabilizing and that the model is converging. This convergence process ceases at the deltaLL value provided below. 

# Typically, the deltaLL would be set as follows:
# Very High Precision and slower time to convergence: 1E-08<->1E-07
# High Precision and slightly faster time to convergence: 1E-06<->1E-05
# Moderate Precision and much faster time to convergence: 1E-04<->1E-03
deltaLL <- 1E-07

# *** The convergence process durint training will cease depending on which param is reached first: deltaLL or numMaxiter. If numMaxiter is set too low, then a 
# model in training may not be able to reach the deltaLL param set, thus its recommended to keep numMaxiter >= 100 to give each model enough iterations to reach
# the deltaLL that was set above. ***

# bepK - background emission probability param for kmer-based models - represents how much additional weight is applied to the background emission probabilities
# of discriminative kmers vs non discriminative kmers. To preferentially score discrim kmer matches higher in models, the background emission probabilities 
# are assigned a biased weighting to discriminative kmers over non-discriminative kmers which are always assigned an initial weighting of 1. 
# Default = 10 (meaning discrim kmers are assigned 10 times the weight to their initial bep compared to the initial bep of non-discrim kmers).
# min value = 1, if set to 1 no biased weighting is applied to the initial background emission probabilities. This does not apply to klen 1 which
# simply uses the background emission probabilities derived from the model.
bepK <- 10

# writeHMM param - should models be written to their own folder in .hmm format (in addition to storing to db)? Represented as a boolean - TRUE/FALSE
# This will allow trained models to be easily exported elsewhere to be used for other sequence dataset classifications however writing models 
# in .HMM will also slightly increase processing time for the kphmm_Train function and will take a great deal more disk space.
writeHMM <- FALSE

# Mixed models param - should mixed kmer length models be trained in addition to individual kmer length models?
# Represented as a boolean - TRUE/FALSE, a mixed model will combine discrim kmers across klengths (excluding klen=1 if it is included for comparison)
# into one model for each group. This param is always disabled when running only one kmer length.
if(length(klenSub)>1){
  mixedModels <- TRUE # Can disable here by setting to FALSE if more than 1 klen but you also dont want mixedModels
} else {
  mixedModels <- FALSE
}

# numSampleTrain param - param to take a small random sample of taxonomic groups and train them - for pipeline testing purposes, set to NA to ignore. 
numSampleTrain <- NA

# Params for kphmm_Train:
# Param 1 - kmerResults - list of discrim kmer results
# Param 2 - kmerData - list of db table names corresponding to each kmer length for frequency and position data
# Param 3 - dfTaxa - dataframe for sequence dataset being used
# Param 4 - klen - numeric for specified kmer length(s)
# Param 5 - numFolds - numeric for specified number of folds used
# Param 6 - numMaxiter - numeric for the maximum number of EM iterations the model must perform before the cycling process is terminated (model setup)
# Param 7 - subclassRank - character for taxonomic ranks to be run
# Param 8 - gapThreshold - numeric for the acceptable proportion of gaps at each position in the model during training
# Param 9 - deltaLL - numeric for the the maximum change in log likelihood between EM iterations before the cycling procedure is terminated (model setup)
# Param 10 - bepK - numeric for background emission probability kmer param, see above for details
# Param 11 - runParallel - boolean indicating whether to enable parallel processing
# Param 12 - numCores - numeric indicating the number of cores to use if runParallel is TRUE, if FALSE, numcores will always equal 1
# Param 13 - futureGlobals - numeric indicating the maximum amount of RAM that is allocated per core when runParallel is TRUE
# Param 14 - doRNGseed - numeric for the seed set for dorng foreach command
# Param 15 - dir - character for the path to write .HMM files to (if enabled)
# Param 16 - dbPath - character for the path to the database being written to
# Param 17 - writeHMM - boolean indicating whether models should also be written to the disk in .hmm format
# Param 18 - mixedModels - boolean indicating whether an additional model that combines all dicrim kmers from all k lengths 
# (excluding klen=1 if it is included) should be trained for each group
# Param 19 - numSampleTrain - numeric for the number of randomly sampled taxonomic groups to train, set to NA to ignore (for pipeline testing purposes)
# Param 20 - runBenchmark - boolean indicating whether to run benchmarking tests

# Model setup and training of PHHM's using all of the specified parameters above
modelData <- kphmm_Train(kmerResults, kmerData, dfTaxa, klen, numFolds, 
                         numMaxiter, subclassRank, gapThreshold, deltaLL, 
                         bepK, runParallel, numCores, futureGlobals, doRNGseed,
                         dir, dbPath, writeHMM, mixedModels, numSampleTrain, 
                         runBenchmark)

# *** Upon completion, trained models will then get written out to a database and will also be written out to files in .hmm format if writeHMM is set to TRUE

# The final output of this function is a list consisting of the following elements per rank/group/fold(if applicable)/klen:
# Group - taxonomic group
# Rank - taxonomic rank
# klen - klength
# fold - fold (if applicable)
# model - db table names of each trained model, this data can be accessed easily, ex: table <- read_dbTable(dbPath, "model_Rhombomys_f2_5")
# success - was model training successful? TRUE/FALSE
# If runBenchmark is TRUE, benchmark data can be accessed easily: table <- read_dbTable(dbPath, "bench_kphmm_Train")

###### 14 - kPHMM Classification (Cross-Validation) ######

# *** The ClassifyCV function will be used when performing cross-validation tests, otherwise the kphmm_ClassifyI function (still in development)  
# should be used when not performing cross-validation tests or when performing classifications independently of training. ***

# Note that the memory overhead is quite high for classifications when using parallel processing, you can adjust numCores higher if your RAM and 
# number of cores is quite high but otherwise you may need to lower the number of cores by 1 or 2

# For testing purposes, to quickly test a random sample of test sequences, set to NA to ignore
numSampleTest <- 50

# Classification strategy options param - must be one of Hierarchical, Multi-Tiered, Post-Processing or MaxScore, default - Hierarchical.
# Must be set to MaxScore if only running at one taxonomic rank.
if(length(subclassRank)==4){
  classStrategy <- "Hierarchical"
} else {
  classStrategy <- "MaxScore"
}

# Must be running at all taxonomic ranks for the following options:

# Hierarchical - sequences are run through models in a purely hierarchical fashion whereby order models are run first, the highest normalized forward score is chosen, family
# level models then get subsetted by the order chosen and only that subset of family models is run. This same process continues for genus and species ranks. 
# This is the fastest approach but also slightly more risky as an incorrect classification made at order level will result in incorrect classifications at all lower levels.

# Multi-Tiered - sequences are run through all Order & Family level models first, a weighted average is then calculated from the normalized forward scores of the order and family models 
# and the viable combination of models (viable meaning a combination of models with a valid taxonomy) with the highest weighted average is chosen. Then following this, genus and 
# species ranks are run hierarchically as described in the Hierarchical method above. This method sacrifices some speed in classifications but has a greater likelihood of making correct 
# classifications at Order and Family level and thus is less risky than Hierarchical.

# Post-Processing - sequences are run through all models at all ranks first and all normalized forward scores are aggregated, a weighted average is then calculated from the normalized forward 
# scores across all levels and the viable combination of models (viable meaning a combination of models with of making classifications since it does not exclude any models from 
# being tested.

# Do not need to be running at all taxonomic ranks for this option:
# MaxScore - the simplest approach, all models are run at each of the rank(s) being run and the model with the highest normalized forward score is chosen. 

# Weight params - weights used in weighted averages for either Multi-Tiered or Post-Processing strategies.
# Default for Multi - c(Order = 0.6, Family = 0.4), default for Post - c(Order = 0.4, Family = 0.3, Genus = 0.2, Species = 0.1).
# Only one set of weights is applied depending on which classification strategy is chosen. 
# Higher ranks are generally weighted higher than lower ranks in these classification strategies but 
# weights can be set as being equal (0.5 for Multi and 0.25 for Post) if a non-weighted average is preferred.
# Set to NA if using Hierarchical or MaxScore.
if(classStrategy == "Multi-Tiered"){
  weights <- c(Order = 0.6, Family = 0.4)
} else if(classStrategy == "Post-Processing"){
  weights <- c(Order = 0.4, Family = 0.3, Genus = 0.2, Species = 0.1)
} else {
  weights <- NA
}

# majorityVoting - an additional boolean param that if TRUE will assign another set of classifications at each rank based on the classifications already
# made across all klengths. If TRUE this voting process will only occur once all classifications for all kmer lengths are performed and can be used with any classification strategy.
# Majority voting means that the classifications across kmer lengths are tallied up based on which groups were assigned, the group with the highest 
# number of votes amongst the various kmer lengths is chosen. Tie breaks are decided by the normalized forward scores assigned, in the event of a tie, the group with highest mean
# normalized forward score is chosen. Default set to TRUE. The classifications made by majority voting will be added alongside the classifications made at individual kmer lengths.
majorityVoting <- TRUE

# Additionally weights are assigned to the voting process depending on kmer length (ex: you want kmer length 1 to be given more weight in voting than other kmer lengths at Species level)
# By default, kmer length 1 results are assigned the highest weight as they are the most reliable currently and already a trusted strategy for PHMM-based classifications in the literature, 
# MixedKlen results (if enabled in training) are assigned the second highest weight since they tend toward higher success rates than other kmer lengths greater than 1 (as of current testing), 
# kmer lengths >= 5 are assigned the third highest weights as they too tend toward higher success rates (as of current testing) and the remaining kmer lengths (2-4) scale downwards in weight 
# by 0.1 increments. Weights can be set anywhere from 0 (no say in voting) to 1 (heaviest weight in voting). Additionally weights are specified per taxonomic rank and thus can be tuned by 
# individual rank. By default all ranks are given the same set of weights but this can be easily modified below. 

# Weights are specified in the data tables below and labeled by their respective ranks:
mvWt <-  list("Order"= data.table("weights" = c(1, 0.9, 0.8, 0.7, 0.6, 0.5),
                                  "klen" = c("1", "MixedKlen", ">=5", "4", "3", "2")),
             "Family"= data.table("weights" = c(1, 0.9, 0.8, 0.7, 0.6, 0.5),
                                  "klen" = c("1", "MixedKlen", ">=5", "4", "3", "2")),
              "Genus"= data.table("weights" = c(1, 0.9, 0.8, 0.7, 0.6, 0.5),
                                  "klen" = c("1", "MixedKlen", ">=5", "4", "3", "2")),
            "Species"= data.table("weights" = c(1, 0.9, 0.8, 0.7, 0.6, 0.5),
                                  "klen" = c("1", "MixedKlen", ">=5", "4", "3", "2")))
mvWt <- rbindlist(mvWt, idcol="rank")
mvWt <- mvWt[(rank %in% subclassRank)]

# Params for kphmm_ClassifyCV:
# Param 1 - kmerData - list of db table names corresponding to each kmer length for frequency and position data
# Param 2 - modelData - list of db table names corresponding to each of the kphmm models trained
# Param 3 - numFolds - numeric for specified number of folds used
# Param 4 - subclassRank - character for taxonomic ranks to be run
# Param 5 - klen - numeric for specified kmer length(s)
# Param 6 - classRank - character for class level taxonomy being run
# Param 7 - dfTaxa - dataframe for sequence dataset being used
# Param 8 - runParallel - boolean indicating whether parallel processing should be used
# Param 9 - numCores - numeric indicating the number of cores to use if runParallel is TRUE, if FALSE, numcores will always equal 1
# Param 10 - futureGlobals - numeric indicating the maximum amount of RAM that can be allocated per core when runParallel is TRUE
# Param 11 - numSampleTest - numeric indicating how many recordIDs (from test dataset) to randomly sample, set to NA to ignore this param
# Param 12 - doRNGseed - numeric for the seed set for dorng foreach command
# Param 13 - dbPath - character for the path to the database being written to
# Param 14 - classStrategy - character for the classification strategy chosen
# Param 15 - weights - numeric for the weights used in weighted average calculations if using Multi-Tiered or Post-Processing options
# Param 16 - runBenchmark - boolean indicating whether to run benchmarking tests
# Param 17 - majorityVoting - boolean indicating whether to use majority voting which will provide an additional set of classifications
# Param 18 - mvWt - data table with each row representing weights assigned for majority voting per rank (if mv is FALSE then this is param is ignored)

# Test each fold combination for each rank on the kPHMM models that have been set up and trained
modelResults <- kphmm_ClassifyCV(kmerData, modelData, numFolds, subclassRank, klen, 
                                 classRank, dfTaxa, runParallel, numCores, futureGlobals, 
                                 numSampleTest, doRNGseed, dbPath, classStrategy, weights, 
                                 runBenchmark, majorityVoting, mvWt)
rm(cData, benchData)

# *** The final output of this function is a list consisting of 3 elements: 
# CorrectIDs - dataframe indicating the correct ids by individual recordID organized by taxonomic rank and klen 
# CorrectID% - dataframe indicating the correct id% by group and fold organized by taxonomic rank and klen
# AccMetrics - dataframe indicating the TP,TN,FP,FN,FPR,TPR,F1,MCC,cohens-kappa & AUC values (calculated from TPR/FPR) organized by taxonomic rank and klen
# ClassifyFail - dataframe indicating instances where classifications for a given sequence failed at a specific rank, klen & fold (if applicable)
# BenchData db table name - if runBenchmark is TRUE, benchmark data can be accessed easily: table <- read_dbTable(dbPath, "bench_kphmm_ClassifyCV") ***

###### 15 - Results Plotting ######

# *** Note that this function only works when cross-validation tests were performed (numFolds >= 3) ***
# *** Also note that the CorrectId% tree plots requires the classifications to have been performed across taxonomic ranks from Order to Species. ***

# The plotResults function is used for plotting results in the form of phylogenetic trees color coded by % correct ids from Order-Species rank, 
# ROC curves (TPR/FPR) per rank color-coded by klen, MCC/F1/cohens-kappa curves (as a function of confidence threshold) per rank, 
# a barplot showing % of groups without discrim kmers per rank/fold/klength, a barplot showing % of groups where model training 
# failed per rank/fold/klength, a barplot showing % of sequences where classifications failed per rank and lastly bar plots of 
# benchmarking data for performance. 

# Plots will be stored in their own folder in the current working directory. 
# This folder is labeled in the format of date_class_Plots, ex: Aug26_Mam_Plots

# dpi setting param for plots written out as png files (default = 100)
dpi <- 100

# Param 1 - kmerData - list generated from the kmer_Find function
# Param 2 - kmerResults - list generated from the kmer_discrimFind function
# Param 3 - modeData - list generated from the kphmm_Train function
# Param 4 - modelResults - list generated from the kphmm_ClassifyCV function
# Param 5 - klen - numeric for the specified kmer lengths
# Param 6 - numFolds - numeric indicating the number of folds being used
# Param 7 - classRank - character for the class rank being run
# Param 8 - subclassRank - character for taxonomic ranks being run
# Param 9 - runBenchmark - boolean indicating whether or not benchmarking was performed
# Param 10 - numCores - numeric for the number of cores that were used in previous functions
# Param 11 - dpi - numeric indicating the desired dpi setting for plots written out as png files (default = 100)

# Run function for plotting results
modelPlots <- plot_ResultsCV(kmerData, kmerResults, modelData, modelResults, 
                             klen, numFolds, classRank, subclassRank, 
                             runBenchmark, numCores, dpi)

# *** The final output of this function is a list consisting of the following elements: 
# Model Plot filepaths (each its own list element) - the file paths of results plots written to disk (in .png format) ***
