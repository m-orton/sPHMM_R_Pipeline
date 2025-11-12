# A strobemer-Based Profile Hidden Markov Model (sPHMM) approach to taxonomic identification of DNA sequence data - Functions

###### Basic Functions ######

# Sequence Import Function to extract out relevant taxonomic data for training and/or cross-validation datasets
fasta_Load <- function(fastaDir, classRank){
  
  # Find files in fastaDir
  fastaFiles <- list.files(path=fastaDir, pattern=paste0("//.fa|.fas",sep=""))
  
  if(length(fastaFiles) == 0){
    stop("No fasta files detected in the directory listed, ensure files are in .fa or .fas format.")
  }
  
  # Load in the file and convert to DNAStringSet (currently only using the first file listed)
  dna <- readDNAStringSet(paste0(fastaDir, "/", fastaFiles, sep=""), format="FASTA")
  
  # Extract names
  dnaNames <- names(dna)
  
  # Define regular expression pattern for name parsing
  pattern <- "[_:;\t ,]"
  
  # Split the input string using strsplit() with the specified pattern above
  dnaNames <- strsplit(dnaNames, pattern)
  
  # Extracting the record id / seq id (this MUST be the first label in the header info of the fasta file)
  dnaNamesId <- sapply(dnaNames, `[`, 1)
  
  # Find which position the class-level taxonomy is in per sequence (since header info may not be consistent across all sequences)
  classPos <- lapply(dnaNames, function(x) which(x == classRank))
  taxaData <- lapply(seq_along(dnaNames), function(i) {
    pos <- classPos[[i]]
    if (length(pos) > 0) {
      start <- pos[1]
      end <- min(start + 5, length(dnaNames[[i]]))
      return(dnaNames[[i]][start:end])
    } else {
      return(NULL)
    }
  })
  
  # For pulling out class level identification and naming the dataframe
  classID <- sapply(taxaData, `[`, 1)
  
  # For order level identification
  ordID <- sapply(taxaData, `[`, 2)
  
  # For family level identification
  famID <- sapply(taxaData, `[`, 3)
  
  # For genus level identification
  genID <- sapply(taxaData, `[`, 4)
  
  # For species level identification
  spID_1 <- sapply(taxaData, `[`, 5)
  spID_2 <- sapply(taxaData, `[`, 6)
  spID <- paste(spID_1, spID_2, sep="_")
  
  # Dataframe with all relevant information for a given class taxonomy in the fasta file
  dfTaxa <- data.frame("recordID" = dnaNamesId, 
                       "dna" = as.character(dna),
                       "classID" = classID,
                       "ordID" = ordID, 
                       "famID" = famID,
                       "genID" = genID,
                       "spID" = spID)
  # rownames ID
  rownames(dfTaxa) <- 1:nrow(dfTaxa)
  
  # Filter out undef labelled records at order/family/genus/species level and na's
  dfTaxa <- dfTaxa %>% filter(ordID != "undef") %>% drop_na()
  dfTaxa <- dfTaxa %>% filter(famID != "undef") %>% drop_na()
  dfTaxa <- dfTaxa %>% filter(genID != "undef") %>% drop_na()
  dfTaxa <- dfTaxa %>% filter(spID != "undef") %>% drop_na()
  
  # Filter out NC recordIDs
  dfTaxa <- dfTaxa %>% filter(recordID != "NC")
  
  # Success message
  cat(green(paste0("Input fasta file successfully loaded for class ", classRank, ".\n")))
  
  # Remove all unnecessary vars created
  rm(dna, dnaNames, dnaNamesId, classPos, taxaData, classID, ordID, famID, genID, spID_1, spID_2, spID)
  return(dfTaxa)
}

# Function for the removal of sequences with gap/n content using the nFilter and gFilter params
# Also filter by a seq length range and minimum number of sequences per species, if any param is NA, ignore that filter
seq_Filter <- function(dfTaxa, nFilter, gFilter, seqLenRange, minNumSpID){
  # N filter
  if(!is.na(nFilter)){
    dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "N|n") / nchar(dfTaxa$dna)) <= nFilter),]
  }
  # Gap filter
  if(!is.na(gFilter)){
    dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "-") / nchar(dfTaxa$dna)) <= gFilter),]
  }
  # Seq length range filtering
  if(!is.na(seqLenRange[1]) & !is.na(seqLenRange[2])){
    dfTaxa <- dfTaxa[which(nchar(dfTaxa$dna) >= seqLenRange[1] & nchar(dfTaxa$dna) <= seqLenRange[2]),]
  }
  # Minimum number of sequences per species
  if(!is.na(minNumSpID)){
    setDT(dfTaxa)
    # Remove species with less than a certain number of sequences (as specified by minNumSpID)
    dfTaxa <- as.data.frame(dfTaxa[, if (.N >= minNumSpID) .SD, by = spID])
  }
}

##### Stratified Sampling Function #####

# Function to create folds based on the number of folds used
stratify_Sampling <- function(dfTaxa, classRank, subclassRank, crossValSeed, numFolds, testNum, trainNum){
  
  # Error for cross-validation params not summing to numFolds if numFolds >= 3
  if(numFolds >= 3 && sum(testNum, trainNum) != numFolds){
    stop("numFolds must be the sum of testNum and trainNum if numFolds is equal to or greater than 3.")
  }
  
  # Class name
  class <- substr(classRank, 1, 3)
  
  # Date for timestamping of new directory
  date <- gsub(' ', '', format(Sys.time(), "%b %d"))
  
  # Top level directory creation for all files with this run
  dir <- paste(date, class, sep="_")
  # dir.create(dir, showWarnings = FALSE)
  
  # Create folds representative of the taxonomy from Order-Species through stratified sampling
  if(numFolds >= 3){
    
    # Update progress  
    cat(green(paste0("Creating fold data for cross validation...\n")))
    
    if(length(subclassRank) == 4){
      
      # a subsetted dataframe for taxonomic rank data which will be used for stratification of sampling
      taxRank <- dfTaxa[c("ordID", "famID", "genID", "spID")]
      
      # Create folds with stratified sampling based on all tax ranks - k here refers to the number of desired strata (or number of k-means clusters)
      # k will end up equal to 4 as long as as each rank has more than one group
      k <- 0
      ifelse(length(unique(taxRank$ordID)) > 1, k <- k + 1, k <- k + 0)
      ifelse(length(unique(taxRank$famID)) > 1, k <- k + 1, k <- k + 0)
      ifelse(length(unique(taxRank$genID)) > 1, k <- k + 1, k <- k + 0)
      ifelse(length(unique(taxRank$spID)) > 1, k <- k + 1, k <- k + 0)
      
      # Stop function if k <= 1
      if(k <= 1){
        stop("Stratified sampling cannot occur using this dataset, not enough unique groups (<= 1 at all levels) are present at each taxonomic rank.")
      }
      
      # Create a stratification variable where k = number of k-means clusters calculated above (split tools package)
      strata <- splitTools::multi_strata(taxRank, k = k)
      
      # list for folds 
      folds <- list()
      # Partition the data using the above stratification variable
      foreach(i=1:numFolds) %do% {
        folds[[i]] <- splitTools::partition(strata, p = c(train = trainNum/numFolds, test = testNum/numFolds), seed = crossValSeed[i], split_into_list = FALSE)
      }
      
      # Create a vector of fold names
      folds <- as.data.frame(folds)
      colnames(folds) <- paste0("Fold", 1:numFolds, sep="")
      
      # merge the fold columns into one column
      folds <- tidyr::unite(folds, col='folds', colnames(folds), sep=';')
      
      # create a taxonomy column for taxonomy data and fold data
      dfTaxa$taxonomy <- paste0(dfTaxa$recordID, ";", dfTaxa$ordID, ";", dfTaxa$famID, ";", dfTaxa$genID, ";", dfTaxa$spID, ";", folds$folds, sep="")
      
    } else {
      
      # Vector for subclassRank column ids
      id <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
      
      # Dataframe isolating only for subclassRank column being run
      taxRank <- dfTaxa[id]
      
      # Create folds with stratified sampling based on one taxrank
      strata <- as.character(taxRank[,1])
      folds <- list()
      foreach(i=1:numFolds) %do% {
        folds[[i]] <- splitTools::partition(strata, p = c(train = trainNum/numFolds, test = testNum/numFolds), seed = seed[i], split_into_list = FALSE)
      }
      
      # Create a vector of fold names
      folds <- as.data.frame(folds)
      colnames(folds) <- paste0("Fold", 1:numFolds, sep="")
      
      # merge the fold columns into one column
      folds <- tidyr::unite(folds, col='folds', colnames(folds), sep=';')
      
      # create a taxonomy column for taxonomy data and fold data
      dfTaxa$taxonomy <- paste0(dfTaxa$recordID, ";", dfTaxa$ordID, ";", dfTaxa$famID, ";", dfTaxa$genID, ";", dfTaxa$spID, ";", folds$folds, sep="")
    }
    
    # update dfTaxa globally
    .GlobalEnv$dfTaxa <- dfTaxa
    
    # add dir var globally
    .GlobalEnv$dir <- dir
    
    # Remove unneeded vars
    rm(folds, strata)
    
  } else if(numFolds == 1){
    # If no cross-val tests to be performed and numFolds is 1
    dfTaxa$taxonomy <- paste0(dfTaxa$recordID, ";", dfTaxa$ordID, ";", dfTaxa$famID, ";", dfTaxa$genID, ";", dfTaxa$spID)
    
    # assign dfTaxa globally
    .GlobalEnv$dfTaxa <- dfTaxa
    
    # add dir var globally
    .GlobalEnv$dir <- dir
    
    # Else if the number of folds is set incorrectly
  } else {
    # Error if numFolds set to a value of 2, less than 1 or a non-integer value
    stop("numFolds must be an integer set to a minimum of 3 to perform cross-validation tests or set to 1 if no cross-validation tests are to be performed.")
  }
}

###### Rcpp Functions ######

# Source and compile all strobemer functions
suppressWarnings(sourceCpp("strobemer_chop.cpp"))
sourceCpp("strobemer_filter.cpp")
sourceCpp("strobemer_filter_unique.cpp")

####### Strobemer Extraction #######

# Function to extract strobemers from sequences at the species level

# find strobes for each sequence of training data with clustering
strobemer_Extract <- function(dtTaxa, nS, kS, wMinS, wMaxS, sType){
  # Sequences of an individual species named with recordID
  seq <- setNames(dtTaxa$dna, dtTaxa$recordID)
  
  # cpp function to find strobemers
  strobeDT <- rbindlist(foreach(i=1:length(seq)) %do% {
    data.table("Strobemer" = chopStrobemer(seq[i], nS, kS, wMinS, wMaxS, sType), "recordID"=names(seq)[i])
  })
  
  # Merging to dtTaxa
  setkey(dtTaxa, recordID)
  setkey(strobeDT, recordID)
  strobeDT <- merge(dtTaxa, strobeDT, allow.cartesian=TRUE)
  strobeDT <- strobeDT[, sID := paste0(recordID, "_", seq_len(.N), sep=""), by = "recordID"]
  
  # Remove unneeded vars
  rm(seq)
  
  # Return strobeDT
  return(strobeDT)
}

####### Stockholm Format Write Function #######

# Function that writes out sequences in stockholm format (as required by HMMER3)
stockholm_Write <- function(strobeDT_Train, fold, dirCyg, classRank){
  # Subset by fold
  strobeClust <- strobeDT_Train[str_detect(names(strobeDT_Train), paste0("f", fold, sep=""))]
  
  # Stockholm format
  strobeClust <- unlist(foreach(i=1:length(strobeClust)) %do% paste0(c("# STOCKHOLM 1.0", 
                                                                       paste0("#=GF ID ", strobeClust[[i]]$clustID[1], sep=""),
                                                                       paste0(strobeClust[[i]]$sID, " ", strobeClust[[i]]$Strobemer),
                                                                       "//"), sep=""))
  
  # Directory for seed files
  dirSeed <- paste0(dirCyg, "strobeClust_f", fold, "_", classRank, ".seed", sep="")
  
  # Writing file out to specified directory
  fileConn <- file(dirSeed)
  writeLines(strobeClust, fileConn)
  close(fileConn)
}

###### Post-Processing of HMMER3 Queries Function ######

# Function to organize and process HMMER3 output data - species level
hmmer3_postProcess <- function(strobeResults, testDict, trainDict){
  # Separate target_name and query name columns into multiple columns
  targetSep <- c("spTrain","fold","sCluster")
  querySep <- c("recordID","ID")
  
  # Create a set of new columns for strobeResults with tstrsplit and filter out hashes
  setDT(strobeResults)
  strobeResults <- strobeResults[clustID != "#", c(targetSep) := tstrsplit(clustID, ";")]
  strobeResults[, c(querySep) := tstrsplit(query, "_")]
  
  # Set key for strobeResults and merge to testDict
  setkey(strobeResults, "recordID")
  strobeResults <- merge(strobeResults, testDict)
  
  # Set key for strobeResults and merge to trainDict
  setkey(strobeResults, "clustID")
  strobeResults <- merge(strobeResults, trainDict, allow.cartesian=TRUE)
  
  # Sort bit score and lowest e-value/bias then select the top ranked sCluster per query
  strobeResults <- strobeResults[order(-score, `E-value`, bias)][,head(.SD, 1), by = c("query")]
  
  # return strobe results
  return(strobeResults)
}

###### Majority Voting Functions #######

# Use a majority voting scheme with weighted sums of the scores to find the most optimal classification for each test recordID, 
# in the event of close scoring between different groups (the param tieThreshold defines what qualifies as a close-score between groups),
# fall back to the raw vote counts to make the classification
classify_MV <- function(dt, group_col, test_col, tieThreshold) {
  
  # Aggregate votes with score weights
  dt_aggregated <- dt[, .(
    weightedScore = sum(score),
    count = .N
  ), by = group_col]
  
  # Find the maximum weighted sum
  maxWeightedScore <- max(dt_aggregated$weightedScore)
  
  # Filter classifications with the maximum weighted sum
  candidates <- dt_aggregated[weightedScore == maxWeightedScore]
  
  # Add back the relevant columns of recordID, fold, and test group
  candidates$recordID <- dt$recordID[1]
  candidates$fold <- dt$fold[1]
  candidates[[test_col]] <- dt[[test_col]][1]
  
  # Fallback to raw vote count in case of ties or close scores
  if (nrow(candidates) > 1) {
    if (diff(range(candidates$weightedScore)) <= tieThreshold) {
      result <- candidates[which.max(count)]
    } else {
      result <- candidates[1]
    }
  } else {
    result <- candidates
  }
  
  # return result df
  return(result)
}

# Function to perform majority vote (with function above) and aggregate majority voting data into one dataframe
aggregate_MV <- function(strobeResults, tieThreshold) {
  # Order
  strobeResults_Ord <- rbindlist(lapply(split(strobeResults, by = c("recordID","fold"), keep.by = TRUE), function(dt) {
    classify_MV(dt, "ordTrain", "ordTest", tieThreshold)
  }))[, rank := "Order"]
  setnames(strobeResults_Ord, c("ordTrain", "ordTest"), c("train","test"))
  
  # Family
  strobeResults_Fam <- rbindlist(lapply(split(strobeResults, by = c("recordID","fold"), keep.by = TRUE), function(dt) {
    classify_MV(dt, "famTrain", "famTest", tieThreshold)
  }))[, rank := "Family"]
  setnames(strobeResults_Fam, c("famTrain", "famTest"), c("train","test"))
  
  # Genus
  strobeResults_Gen <- rbindlist(lapply(split(strobeResults, by = c("recordID","fold"), keep.by = TRUE), function(dt) {
    classify_MV(dt, "genTrain", "genTest", tieThreshold)
  }))[, rank := "Genus"]
  setnames(strobeResults_Gen, c("genTrain", "genTest"), c("train","test"))
  
  # Species
  strobeResults_Sp <- rbindlist(lapply(split(strobeResults, by = c("recordID","fold"), keep.by = TRUE), function(dt) {
    classify_MV(dt, "spTrain", "spTest", tieThreshold)
  }))[, rank := "Species"]
  setnames(strobeResults_Sp, c("spTrain", "spTest"), c("train","test"))
  
  # Aggregate data from all ranks together into one df
  mv_Results <- rbind(strobeResults_Ord, strobeResults_Fam, strobeResults_Gen, strobeResults_Sp)
  
  # Remove unneeded vars
  rm(strobeResults_Ord, strobeResults_Fam, strobeResults_Gen, strobeResults_Sp)
  
  # Return majority vote results
  return(mv_Results)
}

###### Accuracy Benchmarking Functions ######

# Function to calculate accuracy metrics - accuracy, prevalence, sens, spec, f1, cohens kappa, logloss, auc
# MLmetrics is used to determine most metrics, metrica is used to determine balanced accuracy & adjusted f1
calculate_Accuracy <- function(strobeResults){
  # If both negative & positive outcomes exist
  if(nrow(strobeResults[outcome==0])>0 && nrow(strobeResults[outcome==1])>0){
    # use ROCR package and KS statistic to find an optimal cutoff value for other accuracy metrics
    pred <- prediction(strobeResults$weightedScore, strobeResults$outcome)
    perf <- performance(pred, "tpr", "fpr")  # True Positive Rate vs. False Positive Rate
    
    # Compute KS statistic
    ks_values <- perf@y.values[[1]] - perf@x.values[[1]]
    best_ks_index <- which.max(ks_values)  # Index of max separation
    best_cutoff <- pred@cutoffs[[1]][best_ks_index]  # Best threshold
    
    # set of prediction confidence values at best confidence cutoff
    conf_best <- ifelse(strobeResults$weightedScore < best_cutoff, 0, 1)
    
    # MLmetrics for all important metrics
    accDT <- data.table(percent_correct = (nrow(strobeResults[outcome == 1]) / nrow(strobeResults)) * 100,
                        auc = MLmetrics::AUC(y_pred = strobeResults$weightedScore, y_true = strobeResults$outcome),
                        ks_stat = MLmetrics::KS_Stat(y_pred = strobeResults$weightedScore, y_true = strobeResults$outcome),
                        logloss = MLmetrics::LogLoss(y_pred = strobeResults$weightedScore, y_true = strobeResults$outcome),
                        opt_cutoff = best_cutoff,
                        acc = MLmetrics::Accuracy(y_pred = conf_best, y_true = strobeResults$outcome),
                        bal_acc = as.numeric(strobeResults %>% balacc(obs = outcome, pred = conf_best, tidy=TRUE)),
                        error_rate = as.numeric(strobeResults %>% error_rate(obs = outcome, pred = conf_best, tidy=TRUE)),
                        ck = as.numeric(strobeResults %>% khat(obs = outcome, pred = conf_best, tidy=TRUE)),
                        f1 = MLmetrics::F1_Score(y_pred = conf_best, y_true = strobeResults$outcome),
                        adj_f1 = as.numeric(strobeResults %>% agf(obs = outcome, pred = conf_best, tidy=TRUE)),
                        precision = MLmetrics::Precision(y_pred = conf_best, y_true = strobeResults$outcome),
                        recall = MLmetrics::Recall(y_pred = conf_best, y_true = strobeResults$outcome),
                        sens = MLmetrics::Sensitivity(y_pred = conf_best, y_true = strobeResults$outcome),
                        spec = MLmetrics::Specificity(y_pred = conf_best, y_true = strobeResults$outcome))
    
    # Replace NAs with 0
    accDT <- accDT %>% replace(is.na(.), 0)
    
    # Remove unneeded vars
    rm(conf_best, best_cutoff, best_ks_index, ks_values, perf, pred)
  }
  
  # If no negative outcomes exist
  if(nrow(strobeResults[outcome==0])==0){
    accDT <- data.table(percent_correct = (nrow(strobeResults[outcome == 1]) / nrow(strobeResults)) * 100,
                        auc = NA,
                        ks_stat = NA,
                        logloss = NA,
                        opt_cutoff = NA,
                        acc = NA,
                        bal_acc = NA,
                        error_rate = NA,
                        ck = NA,
                        f1 = NA,
                        adj_f1 = NA,
                        precision = NA,
                        recall = NA,
                        sens = NA,
                        spec = NA)
  }
  
  # If no positive outcomes exist
  if(nrow(strobeResults[outcome==1])==0){
    accDT <- data.table(percent_correct = 0,
                        auc = NA,
                        ks_stat = NA,
                        logloss = NA,
                        opt_cutoff = NA,
                        acc = NA,
                        bal_acc = NA,
                        error_rate = NA,
                        ck = NA,
                        f1 = NA,
                        adj_f1 = NA,
                        precision = NA,
                        recall = NA,
                        sens = NA,
                        spec = NA)
  }
  
  # Return dt
  return(accDT)
}

# Main function to determine accuracy for each rank
benchmark_Accuracy <- function(strobeResults){
  # Order
  strobeResults_Ord <- strobeResults[rank == "Order"]
  strobeResults_Ord$outcome <- ifelse(strobeResults_Ord$test == strobeResults_Ord$train, 1, 0)
  
  # Apply platt scaling to raw weighted scoree values
  platt_model <- suppressWarnings(glm(strobeResults_Ord$outcome ~ strobeResults_Ord$weightedScore, family = binomial))
  platt_probs <- predict(platt_model, newdata = data.frame(strobeResults_Ord$weightedScore), type = "response")
  strobeResults_Ord$weightedScore <- platt_probs
  
  # Determine accuracy with calculate_Accuracy function
  results_Ord <- calculate_Accuracy(strobeResults_Ord)
  results_Ord$rank <- "Order"
  
  # Family
  strobeResults_Fam <- strobeResults[rank == "Family"]
  strobeResults_Fam$outcome <- ifelse(strobeResults_Fam$test == strobeResults_Fam$train, 1, 0)
  
  # Apply platt scaling to raw weighted score values
  platt_model <- suppressWarnings(glm(strobeResults_Fam$outcome ~ strobeResults_Fam$weightedScore, family = binomial))
  platt_probs <- predict(platt_model, newdata = data.frame(strobeResults_Fam$weightedScore), type = "response")
  strobeResults_Fam$weightedScore <- platt_probs
  
  # Determine accuracy with calculate_Accuracy function
  results_Fam <- calculate_Accuracy(strobeResults_Fam)
  results_Fam$rank <- "Family"
  
  # Genus
  strobeResults_Gen <- strobeResults[rank == "Genus"]
  strobeResults_Gen$outcome <- ifelse(strobeResults_Gen$test == strobeResults_Gen$train, 1, 0)
  
  # Apply platt scaling to raw weighted score values
  platt_model <- suppressWarnings(glm(strobeResults_Gen$outcome ~ strobeResults_Gen$weightedScore, family = binomial))
  platt_probs <- predict(platt_model, newdata = data.frame(strobeResults_Gen$weightedScore), type = "response")
  strobeResults_Gen$weightedScore <- platt_probs
  
  # Determine accuracy with calculate_Accuracy function
  results_Gen <- calculate_Accuracy(strobeResults_Gen)
  results_Gen$rank <- "Genus"
  
  # Species
  strobeResults_Sp <- strobeResults[rank == "Species"]
  strobeResults_Sp$outcome <- ifelse(strobeResults_Sp$test == strobeResults_Sp$train, 1, 0)
  
  # Apply platt scaling to raw weighted score values
  platt_model <- suppressWarnings(glm(strobeResults_Sp$outcome ~ strobeResults_Sp$weightedScore, family = binomial))
  platt_probs <- predict(platt_model, newdata = data.frame(strobeResults_Sp$weightedScore), type = "response")
  strobeResults_Sp$weightedScore <- platt_probs
  
  # Determine accuracy with calculate_Accuracy function
  results_Sp <- calculate_Accuracy(strobeResults_Sp)
  results_Sp$rank <- "Species"
  
  # Aggregate results from all ranks together into one df
  accResults <- rbind(results_Ord, results_Fam, results_Gen, results_Sp)
  
  # Remove unneeded vars
  rm(strobeResults_Ord, strobeResults_Fam, strobeResults_Gen, strobeResults_Sp,
     results_Ord, results_Fam, results_Gen, results_Sp, platt_probs, platt_model)
  
  # Return results
  return(accResults)
}
