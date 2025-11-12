# *** Ensure R version is updated to at least 4.4.1 prior to using this pipeline ***
# *** Ensure HMMER3 is installed prior to using this pipeline (v 3.3.2 or higher) ***

###### 1 - Working Directory ######

# *** Working directory must include these six files: 
# sPHMM_Pipeline_Train&CV.R, sPHMM_Functions.R, strobemer_extract.cpp, strobemer_chop.cpp, 
# strobemer_extract.h, strobemer_filter.cpp, strobemer_filter_unique.cpp
# and must also include a folder containing the fasta files to be run on this pipeline ex: seqData/
# ex: setwd("C:/sPHMM_testrun) ***
setwd("D/sPHMM_Pipeline_May12") # This path will need to be changed

# Directory to be used for all files used by the HMMER3 analysis
# HMMER3 is currently run on cygwin, thus all files used by HMMER3 are 
# written to the cygwin/home/usr directory however HMMER3 can be run natively on MAC or Linux
# and the directory can just be set to the R working directory above
dirHMMER3 <- "D:/cygwin64/home/Matthew/" # This path will need to be changed

###### 2 - Packages & Functions for Pipeline Functionality ######

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

# List of CRAN packages used in this pipeline
packages <- c("ape", "foreach", "readr", "tidyverse", 
              "plyr", "data.table", "splitTools", "crayon", 
              "MASS", "progressr", "Rcpp", "MLmetrics", 
              "rstatix", "ROCR", "metrica", "future", 
              "doFuture", "gtools", "pryr", "ggplot2", 
              "viridis", "purrr", "ggforce")

# Install (if not installed) and load CRAN packages
for (i in 1:length(packages)) {
  if (!require(packages[i], character.only = TRUE)) {
    install.packages(packages[i], dependencies = TRUE)
    library(packages[i], character.only = TRUE)
  }
}

# One var for all packages used
packages <- sort(append(c("Biostrings"), packages))

# Functions get loaded from a separate R script and must be loaded before the analysis proceeds
source("sPHMM_Functions.R")

###### 3 - Warning/Message Settings ######

# To suppress summarize info in console for the dplyr package
options(dplyr.summarise.inform = FALSE)

# To silence all warnings globally (not recommended since some of the functions below use warnings)
# options(warn=-1)

# Alternatively to turn warnings back on globally if disabled
# options(warn=0)

###### 4 - Global Seed ######

# Global seed value to be used (hmmbuild and nhmmmscan commands for HMMER3)
globSeed <- 42

###### 5 - Benchmarking ######

# ***This section currently does nothing as im working on the code for comparison to standard PHMMs

# if TRUE, will also train standard PHMM models on each group and rank to make comparisons in accuracy
# when compared to the strobemer-based models. 
runBenchmark <- TRUE

##### 6 - Data Import #######

# First name a class to be used for training/testing from the provided fasta file, 
# ex: classRank <- "Arthropoda"

# *** This pipeline currently assumes only one class level taxonomy is being run at one time ***
classRank <- "Branchiopoda"

# *** The sequence data that is imported must be in .fa or .fas format
# and the taxonomic name data labeling in the header info must be in this ordering:
# RecordID/SeqID must be first in the header info, then anywhere after in the header info Class;Order;Family;Genus;Species
# separated by either semi-colon, colon, underscore, space, tab or comma.
# examples of header info compatible with this pipeline with different separator chars: 
# RecordID;...;Class;Order;Family;Genus;Species OR RecordID;Class;Order;Family;Genus;Species OR RecordID;...;Class;Order;Family;Genus;Species;... 
# RecordID,...,Class,Order,Family,Genus,Species OR RecordID,Class,Order,Family,Genus,Species OR RecordID,...,Class,Order,Family,Genus,Species,... ***

# fastaDir - the name of the folder where the fasta file(s) are located within the current working directory ex: seqData/
fastaDir <- "seqData_Branch/"

# Param 1 - fastaDir - name of the folder where the fasta file(s) are located
# Param 2 - classRank - character for the name of the class being run
dfTaxa <- fasta_Load(fastaDir, classRank)

###### 7 - Taxonomic Rank Selection ######

# Var for the taxonomic ranks being run
subclassRank <- c("Order","Family","Genus","Species")

###### 8 - Sequence Filtering ######

# *** Currently only testing pipeline with approx. standard length barcode sequences of 600-700 bp
# but will be expanded to variable length sequences of greater difference in length soon. ***

# Threshold for proportion of N's present in each sequence - Default of 1% (0.01) of total sequence length or less. 
# Set to NA to ignore
nFilter <- 0.01
# Threshold for proportion of gaps present in each sequence - Default of 1% (0.01) of total sequence length or less.
# Set to NA to ignore
gFilter <- 0.01

# Sequence length range to consider for train and test datasets
# Set both values to NA to ignore
seqLenRange <- c(600, 700)

# Minimum number of sequences a species must posses to be used in the analysis 
# (default set to 20 to ensure enough sequences are available at species level for cross-validation)
minNumSpID <- 20

# Param 1 - dfTaxa - loaded sequence dataset 
# Param 2 - nFilter - numeric that represents the proportion of N's present in a sequence, if higher than this threshold, these sequences will be filtered out
# Param 3 - gFilter - numeric that represents the proportion of gap's present in a sequence, if higher than this threshold, these sequences will be filtered out
# Param 4 - seqLenRange -  numeric for the sequence length range to consider for train datasets (and test datasets if performing cross-validation)
# Param 5 - minNumSpID - minimum number of sequences a species must posses to be used in the analysis
dfTaxa <- seq_Filter(dfTaxa, nFilter, gFilter, seqLenRange, minNumSpID)

# ***Insecta only - capping the total number of sequences per species at 50 (randomly sampled from each species)
#setDT(dfTaxa)
#dfTaxa <- dfTaxa[, .SD[sample(.N, min(.N, 50))], by = spID]
#dfTaxa <- data.frame(dfTaxa)
#gc()

###### 9 - Cross-Validation Setup ######

# Creating folds for cross-validation. Subdividing into folds is intended for cross-validation tests

# numFolds - param for total number of folds to create, default set to 5 folds, must be set to a minimum of 3
numFolds <- 5

# testNum, trainNum - params for what proportion of folds to be used for training and testing
# These params must be integers summing to numFolds. These params are only relevant for cross-validation and ignored otherwise
# ex: if numFolds = 5: default proportions of test/train are testNum = 2, trainNum = 3
testNum <- 2
trainNum <- 3

# Set a seed (using global seed) for each fold in the stratified sampling process when performing the cross-validation step in the stratify_Sampling function
crossValSeed <- (globSeed+1:numFolds) - 1

# Stratified Sampling of dfTaxa dataframe - if not performing cross-validation tests, dfTaxa will not be assigned folds

# For stratified sampling, the following params are needed:
# Param 1 - dfTaxa - sequence dataset being used
# Param 2 - classRank - character for the class-level taxonomy currently being run
# Param 3 - subclassRank - character for the all subclass ranks being run
# Param 4 - crossValSeed - numeric vector for the seed set to each fold in cross-validation
# Params 5-7 - numeric vectors for cross-validation params - numFolds, testNum & trainNum
stratify_Sampling(dfTaxa, classRank, subclassRank, crossValSeed, 
                  numFolds, testNum, trainNum)

# Garbage collection
gc()

###### 10 - Parallel Processing and Progress Bar ######

# Set to TRUE to enable parallel processing for model setup/training and testing/cross-validation.
# *** Recommended for faster runtimes but comes at the cost of higher overall memory usage. ***
runParallel <- TRUE

# if TRUE, set global settings for parallel processing
if(runParallel == TRUE){
  # Amount of RAM to allocate to each core for parallel processing (1.5 Gb default)
  futureGlobals <- 1500*1024^2
  # numCores param - Default set to available cores if runParallel is TRUE
  # If set to FALSE, numCores will always be set to 1. Recommended to subtract 1 core if 
  # runParallel is TRUE and this program is being run locally
  numCores <- availableCores() #- 1
  options(future.globals.maxSize = futureGlobals)
  registerDoFuture()
  plan(multisession, workers=numCores)
  # False, numCores is always set to 1
} else {
  plan(sequential)
  numCores <- 1
}

# Initialize progress bar
handlers("progress")

###### 11 - Important Params ######

# Params for extraction of strobemers (both training sequences and testing sequences)
nS <- 8 # number of strobes
kS <- 20 # kmer size to use
wMinS <- 21 # Min window size for strobes
wMaxS <- 75 # Max window size for strobes
sType <- 1 # 1 for ministrobes, 2 for randstrobes, 3 for hybridstrobes. Using ministrobes by default, randstrobes and hybridstrobes are not recommended

# Total strobemer length calculation
strobeLength <- kS * nS

# Params for strobemer filtering & clustering

# Maximum number of top ranked strobemers per species and fold
# Higher values of maxS will mean longer training and classification times but usually greater accuracy,
# 50 seems to be a good middle ground to start with, going lower than 50 means there may be an inadequate number of strobemer clusters
# to train models on
maxS <- 50
# minimum number of sequences that must be present in a cluster for a model to be trained on that cluster (cannot be set lower than 2)
minClustSize <- 2
# cutoff threshold for clustering of strobemer sequences (excluding the minimizer region), represents the maximum number
# of mismatches allowed between sequences within a cluster, default set as being 15
tC <- 15

# Params for model testing (test sequence sampling, strobemer sampling and majority voting)

# Param for how many records to sample per fold and species (test records are stratified across species in sampling to ensure an even distribution)
recordSampSize <- 3
# Strobemer sampling interval across all strobemers extracted for a given test sequence, default set at 5
# for example, if set at 10, sample 1 strobemer every 10 strobemers across all strobemers extracted from a given test sequence
strobeSampInterval <- 5
# a threshold that identifies how close scores have to be in majority voting before a fallback to raw vote counts is implemented
tieThreshold <- 0.1

# HMMER3 reporting bit score threshold, default set at 20, must be set carefully,
# if set too low, HMMER3 will report and write-out hundreds to thousands of hits slowing down its program runtime and increasing disk usage, 
# if set too high, HMMER3 may not be able to report any results for certain test sequences
bitThresholdReport <- 20

###### 12 - Strobemer Extraction ######

# Set up dtTaxaSp dataframe for finding strobemers
setDT(dfTaxa)

# Create naming vectors for the taxonomy column of dfTaxaSp
ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
nameVecFolds <- paste0("Fold", 1:numFolds, sep="")
nameVecTax <- c("recordID", ids, nameVecFolds)

# Create a set of fold columns for dfTaxa with separate
dfTaxa <- dfTaxa %>%
  ungroup() %>%
  separate(taxonomy, nameVecTax, ";")

# Grab the necessary fold data first before proceeding
foldCols <- paste0("Fold", 1:numFolds, sep="")

# Training sequences strobemer extraction

# Train data
dtTaxa_Train <- rbindlist(foreach(i=1:numFolds) %do% {
  # Subset the data into different folds
  dfTaxa[which(dfTaxa[,i+7] == "train"),] %>%
    mutate(folds = paste0("f", i))
}) %>% select(!all_of(foldCols)) %>% as.data.table()

# Use the find strobemers function to extract from each fold
strobeDT_Train <- with_progress({
  p <- progressor(along = 1:numFolds, offset = 1)
  rbindlist(lapply(split(dtTaxa_Train, by = c("folds"), keep.by = TRUE), function(dt) {
    p("Extracting strobemers...")
    strobemer_Extract(dt, nS, kS, wMinS, wMaxS, sType)
  }))
})

# Garbage collection
gc()

# Restrict columns to only the most relevant ones
strobeDT_Train <- strobeDT_Train[,c("recordID","spID","genID","famID","ordID","Strobemer","folds")]

# Garbage collection
gc()

# Test sequence strobemer extraction and sampling

# Test data
dtTaxa_Test <- rbindlist(foreach(i=1:numFolds) %do% {
  # Subset the data into different folds
  dfTaxa[which(dfTaxa[,i+7] == "test"),] %>%
    mutate(folds = paste0("f", i))
}) %>% select(!all_of(foldCols)) %>% as.data.table()

# Split by spID and folds
strobeDT_Test <- split(dtTaxa_Test, by = c("spID","folds"), keep.by = TRUE)
# Set seed for test record sampling
set.seed(globSeed)
# sample records (number set by recordSampleSize) per spID and fold
strobeDT_Test <- rbindlist(foreach(i=1:length(strobeDT_Test)) %do% strobeDT_Test[[i]][sample(seq_along(strobeDT_Test[[i]]), size = recordSampSize)])
# Remove NA rows
strobeDT_Test <- strobeDT_Test[!is.na(strobeDT_Test$dna), ]
# Split by recordID 
strobeDT_Test <- split(strobeDT_Test, by = "recordID", keep.by = TRUE)

# Extract strobemers from test sequences sampled
strobeDT_Test <- with_progress({
  p <- progressor(along = 1:length(strobeDT_Test), offset = 1)
  rbindlist(foreach(i=1:length(strobeDT_Test)) %do% {
    p("Extracting test strobemers...")
    uniqFolds <- unique(strobeDT_Test[[i]]$folds)
    rbindlist(foreach(j=1:length(uniqFolds)) %do% {
      strobemer_Extract(strobeDT_Test[[i]][(folds == uniqFolds[j])], nS, kS, wMinS, wMaxS, sType)
    })
  })
})

# Sample strobemers from test sequences (sampling determined by strobeSampInterval param)
strobeDT_Test <- split(strobeDT_Test, by = c("recordID","folds"), keep.by = TRUE)
strobeDT_Test <- with_progress({
  p <- progressor(along = 1:length(strobeDT_Test), offset = 1)
  rbindlist(foreach(i=1:length(strobeDT_Test)) %do% {
    p("Sampling test strobemers...")
    strobeDT_Test[[i]][, .SD[seq(1, .N, by = strobeSampInterval)]]
  })
})

# Remove unneeded vars and garbage collection
rm(uniqFolds, foldCols, ids, nameVecTax, nameVecFolds)
gc()

###### 13 - Strobemer Filtering & Clustering ######

# Finding all unique strobemers per species (and then splitting by fold) with the rcpp function findUniqueSequences, this ensures all strobemers
# present are unique to a particular species and are not shared between species, strobemers are particularly good at discriminating at species level
strobeDT_Train <- with_progress({
  p <- progressor(along = 1:numFolds, offset = 1)
  rbindlist(lapply(split(strobeDT_Train, by = c("folds"), keep.by = TRUE), function(dt) {
    p("Finding unique strobemers...")
    strobemer_filterUniqueSp(dt)
  }), idcol="folds")
})

# Garbage collection
gc()

# Split by fold and then split again by spID
strobeDT_Train <- split(strobeDT_Train, by = c("folds"), keep.by = TRUE)
strobeDT_Train <- lapply(strobeDT_Train, function(subset) {
  split(subset, subset$spID)
})

# Rank strobemers by using a metric defined here as sCoverage (strobemer Coverage) - the proportion of recordIDs/seqIDs within a species
# that posses a specific strobemer. Strobemers that posses a high sCoverage are shared by a larger proportion of the seqID/recordIDs 
# within a species (and fold) and thus are more likely to be more useful in classification and clustering. This uses the rcpp function strobemer_Filter.
strobeDT_Train <- with_progress({
  p <- progressor(along = 1:numFolds, offset = 1)
  rbindlist(foreach(i=1:numFolds) %do% {
    p("Sorting & filtering strobemers...")
    rbindlist(foreach(j=1:length(strobeDT_Train[[i]])) %do% {
      foldName <- as.character(strobeDT_Train[[i]][[j]]$folds[1])
      spIDName <- as.character(strobeDT_Train[[i]][[j]]$spID[1])
      strobemer_FilterSp(strobeDT_Train[[i]][[j]], foldName, spIDName, maxS)
    })
  })
})

# Garbage collection
rm(spIDName, foldName)
gc()

# Taxonomy dictionary to reference taxonomy for training dataset
taxDict <- dfTaxa[,c("ordID","famID","genID","spID")]
setDT(taxDict)
taxDict <- unique(taxDict)
setkey(taxDict, spID)

# Merge to strobeDT dt
setkey(strobeDT_Train, spID)
strobeDT_Train <- merge(strobeDT_Train, taxDict)

# Extract the downstream sub sequence of the strobemer (defined as the sequence kS+1 bp downstream of the full strobemer sequence)
strobeDT_Train$strobeSubseq <- substr(strobeDT_Train$Strobemer, kS+1, nchar(strobeDT_Train$Strobemer))

# Create a unique identifier for each strobemer called sID
strobeDT_Train <- strobeDT_Train[, sID := paste0(spID, ";", folds, ";", seq_len(.N), sep=""), by = c("spID", "folds")]

# Split by fold and then split again by spID
strobeDT_Train <- split(strobeDT_Train, by = c("folds"), keep.by = TRUE)
strobeDT_Train <- lapply(strobeDT_Train, function(subset) {
  split(subset, subset$spID)
})

# Clustering strobemers with the stringDist and cuttree functions
strobeDT_Train <- with_progress({
  p <- progressor(along = 1:numFolds, offset = 1)
  rbindlist(foreach(i=1:numFolds) %dofuture% {
    p("Clustering strobemers...")
    rbindlist(foreach(j=1:length(strobeDT_Train[[i]])) %do% {
      if(length(strobeDT_Train[[i]][[j]]$Strobemer)>=2){
        # Calculate pairwise Hamming distances (hamming is chosen because all strobemers will be of uniform length)
        dist_matrix <- pwalign::stringDist(DNAStringSet(setNames(strobeDT_Train[[i]][[j]]$strobeSubseq, strobeDT_Train[[i]][[j]]$sID)), method = "hamming")
        
        # cluster sequences using cutree and hclust functions with specified clustering threshold
        sClusters <- split(setNames(strobeDT_Train[[i]][[j]]$strobeSubseq, strobeDT_Train[[i]][[j]]$sID), 
                           cutree(hclust(as.dist(dist_matrix), method = "average"), k = NULL, h = tC))
        
        # Convert to data.table with a unique cluster ID for each cluster, sID and the strobeSubseqs that comprise each cluster
        sClusters <- rbindlist(
          lapply(names(sClusters), function(list_name) {
            data.table(
              clustID = list_name,
              sID = names(sClusters[[list_name]]),
              strobeSubseq = unname(sClusters[[list_name]])
            )
          }),
          use.names = TRUE
        )
        # Finally merge back to strobeDT_Train and modify the clustID with species id and fold
        sClusters <- merge(strobeDT_Train[[i]][[j]], sClusters, by = c("sID","strobeSubseq"))[, clustID := paste0(spID, ";", folds, ";", clustID, sep="")]
      } else {
        NULL
      }
    })
  })
})

# Garbage collection
gc()

# Remove clusters with less than a specified number of sequences (minClustSize)
strobeDT_Train <- strobeDT_Train[, if (.N >= minClustSize) .SD, by = clustID]

# Split by unique cluster
strobeDT_Train <- split(strobeDT_Train, by = c("clustID"), keep.by = TRUE)

# Garbage collection
rm(taxDict)
gc()

###### 14 - Model Training Setup (HMMER3) ######

# Lists for hmm commands
hmmComm <- list()
dirSh <- list()

# Write out strobemers per unique fold
foreach(i=1:numFolds) %do% stockholm_Write(strobeDT_Train, i, dirHMMER3, classRank)

# HMMER3 commands used to train phmm's and build a phmm database
hmmBuild <- foreach(i=1:numFolds) %do% paste0("hmmbuild --cpu ", numCores, " --seed ", globSeed, " strobeClust_f",
                                              i, "_", classRank, ".hmm", " strobeClust_f", i, "_", classRank, ".seed", sep="")
hmmPress <- foreach(i=1:numFolds) %do% paste0("hmmpress strobeClust_f", i, "_", classRank, ".hmm", sep="")
hmmComm[[1]] <- append(hmmBuild, hmmPress)

# Shebang line for .sh file is appended
hmmComm[[1]] <- append("#!/bin/bash", hmmComm[[1]])

# Filepath for sh file containing all commands to run
dirSh <- list(paste0(dirHMMER3, "run_HMM_Train", "_", classRank, ".sh", sep=""),
              paste0(dirHMMER3, "run_HMM_Test", "_", classRank, ".sh", sep=""))

# Write out the sh script containing all training and testing commands for HMMER3
con <- file(dirSh[[1]])
if( !isOpen(con = con, rw = "wb") ) { open( con, open = "wb" ) }
writeLines(unlist(hmmComm[[1]]), con)
close(con)

# Remove uneeded vars
rm(hmmPress, hmmBuild)

# command to execute the .sh script that contains all HMMER3 commands
writeClipboard(dirSh[[1]])

###### 15 - Model Querying Setup (HMMER3) ######

# Write out to fasta file for testing with HMMER3 (by fold)
foreach(i=1:numFolds) %do% {
  writeXStringSet(DNAStringSet(setNames(strobeDT_Test[folds == paste0("f", i, sep="")]$Strobemer, 
                                        strobeDT_Test[folds == paste0("f", i, sep="")]$sID)), 
                  paste0(dirHMMER3, "strobeTest_f", i, "_", classRank, sep=""))
}

# HMMER3 commands to run for nhmmscan
hmmScan <- foreach(i=1:numFolds) %do% paste0("nhmmscan --cpu ", numCores, 
                                             " -T ", bitThresholdReport,
                                             " --seed ", globSeed, 
                                             " --noali --dfamtblout strobeResults_f", i,  
                                             "_", classRank, ".txt", " strobeClust_f", i, "_", classRank, ".hmm", 
                                             " strobeTest_f", i, "_", classRank, sep="")

# Combining commands for .sh scripts
hmmComm[[2]] <- append("#!/bin/bash", hmmScan)

# Remove uneeded vars
rm(hmmScan)

# Write out the sh script containing all training and testing commands for HMMER3
con <- file(dirSh[[2]])
if( !isOpen(con = con, rw = "wb") ) { open( con, open = "wb" ) }
writeLines(unlist(hmmComm[[2]]), con)
close(con)

# command to execute the .sh script that contains all HMMER3 commands
writeClipboard(dirSh[[2]])

####### 16 - Post-Processing of Model Queries & Classification #######

# Create a dictionary of test records along with their taxonomic IDs 
testDict <- strobeDT_Test[,c("recordID", "ordID","famID","genID","spID")]
colnames(testDict)[2:5] <- c("ordTest","famTest","genTest","spTest")
setDT(testDict)
testDict <- unique(testDict)
setkey(testDict, "recordID")

# Create a dictionary of train records along with their taxonomic IDs 
trainDict <- rbindlist(strobeDT_Train)[,c("clustID","ordID","famID","genID")]
colnames(trainDict)[2:4] <- c("ordTrain","famTrain","genTrain")
setDT(trainDict)
trainDict <- unique(trainDict)
setkey(trainDict, "clustID")

# Read in results from HMMER3 by fold
strobeResults <- with_progress({
  p <- progressor(along = 1:numFolds, offset = 1)
  rbindlist(foreach(i=1:numFolds) %dofuture% {
    p("Processing HMMER3 results...")
    
    # Read in HMMER3 results
    strobeResults <- read.table(paste0(dirHMMER3, "strobeResults_f", i, "_", classRank, ".txt", sep=""), quote="\"")[,c(1,3,4,5,6)]
    colnames(strobeResults) <- c("clustID","query","score","E-value","bias")
    
    # Use function hmmer3_postProcess function to organize strobe results data
    strobeResults <- hmmer3_postProcess(strobeResults, testDict, trainDict)
  })
})

# Garbage collection
gc()

# Apply the majority vote function at each taxonomic rank and aggregate results across taxonomic ranks
strobeResults_MV <- aggregate_MV(strobeResults, tieThreshold)

###### 17 - Accuracy Benchmarking ######

# Benchmarking of important accuracy metrics with benchmark accuracy function (will benchmark individually for each taxonomic rank)
accResults <- benchmark_Accuracy(strobeResults_MV)

# Add all important params to results to determine the effect they might have
# classRank, minclustSize, maxS, strobeSampInterval, recordSampSize, ns, ks, wMinS, wMaxS, sType, minCntMV, seqLenRange
accResults$classID <- classRank
accResults$seqLenMin <- round(min(nchar(dfTaxa$dna)), 0)
accResults$seqLenMed <- round(median(nchar(dfTaxa$dna)), 0)
accResults$seqLenMax <- round(max(nchar(dfTaxa$dna)), 0)
accResults$seqLenRange <- round(max(nchar(dfTaxa$dna)), 0) - round(min(nchar(dfTaxa$dna)), 0)
accResults$recSampSize <- recordSampSize
accResults$totalNumRecords <- length(unique(strobeDT_Test$recordID))
accResults$strobeSampInterval <- strobeSampInterval
accResults$minClustSize <- minClustSize
accResults$tC <- tC
accResults$maxS <- maxS
accResults$strobeLength <- strobeLength
accResults$numStrobes <- nS
accResults$kSize <- kS
accResults$minWindow <- wMinS
accResults$maxWindow <- wMaxS
accResults$strobeType <- "ministrobe"

# Create another directory for the results
dir.create("sPHMM_Results", showWarnings = FALSE)

# Output results to tsv format
write_tsv(accResults, paste0("sPHMM_Results/accResults_", classRank, "_", recordSampSize, "_",
                             (seqLenRange[2] - seqLenRange[1]), "_", strobeSampInterval, "_", 
                             minClustSize, "_", tC, "_", maxS, "_", nS, "_", kS, "_", 
                             wMinS, "_", wMaxS, "_ministrobes.tsv", sep=""))

###### 18 - Results Plotting ######

# Plot a heatmap for all accuracy metrics across different taxonomic groups

# First read in tsv result files
tsv_results <- list.files(path = "D:/sPHMM_Pipeline/sPHMM_Results", pattern = "\\.tsv$", full.names = TRUE) # this path needs to be changed as well

# Read and combine them into one data frame
heat_results <- tsv_results %>%
  map_dfr(~ read_tsv(.x, col_types = cols()))

# Select desired columns
heat_results <- heat_results[,c(1:4,6:7,9:17,23)]

# Add number records tested to classID names
heat_results$classID <- paste0(heat_results$classID, " (nTest=", heat_results$totalNumRecords, ")", sep="")
heat_results <- heat_results[,c(1:15)]

# Convert percent correct back to decimal and represent as fraction correct
heat_results <- heat_results %>% 
  mutate(fraction_correct = percent_correct / 100)
heat_results <- heat_results[,c(16,2:15)]

# Convert to long format and convert rank and metric to factors (for proper ordering in plot)
heat_results <- heat_results %>%
  pivot_longer(cols = c(fraction_correct:spec), # Select numeric columns
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(rank = factor(rank, levels = c("Species","Genus","Family","Order"))) %>%
  mutate(classID = gsub("Actinopteri", "Actinopterygii", classID)) %>% # Only relevant if running fishes
  mutate(Metric = factor(Metric, levels = c("fraction_correct", "auc", "ks_stat",  
                                            "opt_cutoff", "acc", "bal_acc", "ck", 
                                            "f1", "adj_f1", "precision", "recall", 
                                            "sens", "spec", "logloss")))

# Factor names for classID
factor_names <- unique(heat_results$classID)

# Convert classID to factor
heat_results$classID <- factor(as.character(heat_results$classID), labels=factor_names)

# New labels for plot
new_labels <- c(
  "fraction_correct" = "Fraction Correct",
  "auc" = "AUC",
  "ks_stat" = "KS Statistic",
  "logloss" = "Log Loss",
  "acc" = "Accuracy",
  "bal_acc" = "Balanced Acc.",
  "ck" = "Cohen's Kappa",
  "f1" = "F1",
  "adj_f1" = "Adjusted F1",
  "precision" = "Precision",
  "recall" = "Recall",
  "sens" = "Sensitivity",
  "spec" = "Specificity",
  "rank" = "Model Rank"
)

# Plot heatmaps for ranks for Order -> Species
heatmap <- list()
foreach(i=1:ceiling(length(factor_names)/8)) %do% {
  heatmap[[i]] <- ggplot(heat_results, aes(x = Metric, y = rank, fill = Value)) +
    geom_tile(color = "black") + # set black border around tiles
    scale_fill_viridis_c(limits = c(0, 1)) +  # Set fixed color scale from 0 to 1
    scale_x_discrete(labels = new_labels) +  # Change x-axis labels
    facet_wrap_paginate(~classID, ncol = 1, nrow = 8, page=i) + # Facet wrap by classID
    theme_dark(base_size = 16) + # Dark theme and text size
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
    labs(title = "", x = "", y = "")
  # Save plot
  ggsave(paste0("sPHMM_Results/acc_heatmap_", i, ".png", sep=""), plot = heatmap[[i]], width = 9, height = 11, dpi = 300,  bg = "white")
}
