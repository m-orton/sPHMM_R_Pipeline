# A kmer-Based Profile Hidden Markov Model (kPHMM) approach to taxonomic identification of DNA sequence data - Functions

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
  
  # Filter out undef labelled records at order level and na's
  dfTaxa <- dfTaxa %>% filter(ordID != "undef") %>% drop_na()
  
  # Filter out NC recordIDs
  dfTaxa <- dfTaxa %>% filter(recordID != "NC")
  
  # Success message
  cat(green(paste0("Input fasta file successfully loaded for class ", classRank, ".\n")))
  
  # Remove all unecessary vars created
  rm(dna, dnaNames, dnaNamesId, classPos, taxaData, classID, ordID, famID, genID, spID_1, spID_2, spID)
  return(dfTaxa)
}

# Function that creates a character vector based on which taxonomic ranks are specified by the user
subclass_Rank <- function(subclassRankOpt, classRank){
  
  # If all ranks are being analyzed
  if(subclassRankOpt == "ALL" | subclassRankOpt == "All"){
    # subclassRank vector for the four taxonomic groups to use:
    subclassRank <- c("Order", "Family", "Genus", "Species")
    cat(green(paste0("Analyses will be conducted at the following taxonomic ranks: Order, Family, Genus and Species for Class ", classRank, ".\n")))
    # If analyses are only being performed at one rank
    # Species
  } else if (subclassRankOpt == "Species" | subclassRankOpt == "species"){
    subclassRank <- c("Species")
    cat(green(paste0("Analyses will be conducted at the following taxonomic rank: ", subclassRankOpt, " for Class ", classRank, ".\n")))
    # Genus
  } else if (subclassRankOpt == "Genus" | subclassRankOpt == "genus"){
    subclassRank <- c("Genus")
    cat(green(paste0("Analyses will be conducted at the following taxonomic rank: ", subclassRankOpt, " for Class ", classRank, ".\n")))
    # Family
  } else if (subclassRankOpt == "Family" | subclassRankOpt == "family"){
    subclassRank <- c("Family")
    cat(green(paste0("Analyses will be conducted at the following taxonomic rank: ", subclassRankOpt, " for Class ", classRank, ".\n")))
    # Order
  } else if (subclassRankOpt == "Order" | subclassRankOpt == "order"){
    subclassRank <- c("Order")
    cat(green(paste0("Analyses will be conducted at the following taxonomic rank: ", subclassRankOpt, " for Class ", classRank, ".\n")))
    # If incorrect argument is specified in the function
  } else {
    stop(paste0("The taxonomic rank specified: ", subclassRankOpt, " does not belong to one of either All, Order, Family, Genus or Species.\n"))
  }
  # Return the subclassRank vector
  return(subclassRank)
}

# Function for the removal of sequences with gap/n content using the nFilter and gFilter params
# Also filter by a seq length range, if any param is NA, ignore that filter
seq_Filter <- function(dfTaxa, nFilter, gFilter, seqLenRange){
  if(!is.na(nFilter)){
    dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "N|n") / nchar(dfTaxa$dna)) <= nFilter),]
  }
  if(!is.na(gFilter)){
    dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "-") / nchar(dfTaxa$dna)) <= gFilter),]
  }
  if(!is.na(seqLenRange[1]) & !is.na(seqLenRange[2])){
    dfTaxa <- dfTaxa[which(nchar(dfTaxa$dna) >= seqLenRange[1] & nchar(dfTaxa$dna) <= seqLenRange[2]),]
  }
}

# Function to sample with the minimum number and maximum number of sequences at each taxonomic group and rank
subclass_Filter <- function(dfTaxa, subclassRank, groupMinNum, groupMaxNum){
  
  # Stop function if groupMinNum is set to less than 2
  if(groupMinNum < 2){
    stop("The value for groupMinNum cannot be less than 2.")
  }
  
  # If all ranks are being analyzed, group by all taxonomic ranks (which really means grouping and filtering by species level) and filter 
  # by groupMinNum & groupMaxNum
  if(length(subclassRank == 4)){
    if(!is.na(groupMaxNum)){
      dfTaxa %>% 
        group_by(ordID, famID, genID, spID) %>% 
        filter(n() >= groupMinNum) %>% 
        dplyr::slice(if(n() > groupMaxNum) sample(1:n(), groupMaxNum) else 1:n())
    } else {
      dfTaxa %>% 
        group_by(ordID, famID, genID, spID) %>% 
        filter(n() >= groupMinNum)
    }
    # If only one rank is being analyzed
  } else {
    # Order
    if(subclassRank == "Order"){
      if(!is.na(groupMaxNum)){
        dfTaxa %>% 
          group_by(ordID) %>% 
          filter(n() >= groupMinNum) %>% 
          dplyr::slice(if(n() > groupMaxNum) sample(1:n(), groupMaxNum) else 1:n())
      } else {
        dfTaxa %>% 
          group_by(ordID) %>% 
          filter(n() >= groupMinNum)
      }
    }
    # Family
    if(subclassRank == "Family"){
      if(!is.na(groupMaxNum)){
        dfTaxa %>% 
          group_by(famID) %>% 
          filter(n() >= groupMinNum) %>% 
          dplyr::slice(if(n() > groupMaxNum) sample(1:n(), groupMaxNum) else 1:n())
      } else {
        dfTaxa %>% 
          group_by(famID) %>% 
          filter(n() >= groupMinNum)
      }
    }
    # Genus
    if(subclassRank == "Genus"){
      if(!is.na(groupMaxNum)){
        dfTaxa %>% 
          group_by(genID) %>% 
          filter(n() >= groupMinNum) %>% 
          dplyr::slice(if(n() > groupMaxNum) sample(1:n(), groupMaxNum) else 1:n())
      } else {
        dfTaxa %>% 
          group_by(genID) %>% 
          filter(n() >= groupMinNum)
      }
    }
    # Species
    if(subclassRank == "Species"){
      if(!is.na(groupMaxNum)){
        dfTaxa %>% 
          group_by(spID) %>% 
          filter(n() >= groupMinNum) %>% 
          dplyr::slice(if(n() > groupMaxNum) sample(1:n(), groupMaxNum) else 1:n())
      } else {
        dfTaxa %>% 
          group_by(spID) %>% 
          filter(n() >= groupMinNum)
      }
    }
  }
}

# Internal function to write out serialized data to a database
.write_dbTable <- function(db_path, table_name, data) {
  # Serialize and compress data
  data <- qserialize(data, preset = "high")
  raw_data <- data.frame("1" = I(list(data)))
  colnames(raw_data)[1] <- table_name
  # Open connection and write serialized data out to db
  con <- dbConnect(RSQLite::SQLite(), db_path)
  on.exit(dbDisconnect(con))
  dbWriteTable(con, table_name, raw_data, row.names = FALSE, overwrite = TRUE)
}

# Internal function to write out serialized data to a database
read_dbTable <- function(dbPath, data) {
  # Open connection and read in serialized data from db
  con <- dbConnect(RSQLite::SQLite(), dbPath)
  on.exit(dbDisconnect(con))
  data <- dbReadTable(con, data)
  # Deserialize data
  qdeserialize(unlist(data[, 1]), use_alt_rep = FALSE, strict = FALSE)
}

# Internal function for scaling values 0-1
.scaleValues <- function(x){(x-min(x))/(max(x)-min(x))}

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

##### Discriminative Kmer Functions #####

# Internal function to load the C++ kmer search code for grabbing the positional data of the kmers (using the Rcpp R package)
# Source for the C++ code
sourceCpp("seqstokmers.cpp")

# Main function to find both the counts (across sequence window specified by the user) and the positional data (cpp script)
# The main function used for kmer extraction function
kmer_Find <- function(dfTaxa, subclassRank, dir, klen, runBenchmark, runParallel, numCores, futureGlobals, doRNGseed){
  
  # Initialize parallel processing settings
  if(runParallel == TRUE){
    options(future.globals.maxSize = futureGlobals)
    registerDoFuture()
    plan(multisession, workers=numCores)
    handlers("progress")
  } else {
    plan(sequential)
    handlers("progress")
  }

  # Path to the SQLite database
  dbPath <- paste0(dir, ".sqlite", sep="")
  
  # Set as global var
  .GlobalEnv$dbPath <- dbPath
  
  # Convert to DNAbin format
  dnaBin <- as.DNAbin(DNAStringSet(dfTaxa$dna))
  names(dnaBin) <- dfTaxa$taxonomy
  
  # klen var without 1
  klen1Remove <- klen[klen!=1]
  
  # Create naming vectors for the taxonomy column of each kData dataframe depending on numFolds & subclassRank
  ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
  if(numFolds >= 3){
    nameVecFolds <- paste0("Fold", 1:numFolds, sep="")
    nameVecTax <- c("recordID", ids, nameVecFolds)
  } else {
    nameVecTax <- c("recordID", ids)
  }
  
  # Creating list of table names for kmer count data
  kmerFiles_c <- foreach(i=1:length(klen1Remove)) %do% paste("raw_klen_", klen1Remove[i], "_cnt", sep="")
  
  # Garbage collection
  gc()
  
  # Perform kmer counts with kcount function from kmer R package for all klengths greater than 1
  # and write to tsv files
  bench_cnt <- with_progress({
    p <- progressor(along = 1:length(klen1Remove), offset = 1)
    foreach(i=1:length(klen1Remove), .options.RNG=doRNGseed, .packages = "RSQLite") %dorng% {
      # Update message
      p(paste0("Performing kmer counts and writing to db", " (klen=", klen1Remove[i], ")", sep=""))
      
      # Bench st
      benchS <- .bench_Func(runBenchmark, klen1Remove[i], "kmer_Find", "kCnt", 0, "Start")

      # use kmer package to grab kmer data for all sequences
      count_matrix <- kmer::kcount(dnaBin, k = klen1Remove[i], residues = "DNA")
      
      # Sum the counts across columns for with Rfast package
      col_totals <- colSums(count_matrix)
      
      # Calculate frequencies for columns
      count_matrix <- count_matrix / col_totals
      
      # Convert to dataframe format
      count_matrix <- as.data.frame(count_matrix) %>% 
        rownames_to_column(var="taxonomy") %>%
        ungroup() %>%
        separate(taxonomy, nameVecTax, sep = ";")
      
      # Write count data to db
      .write_dbTable(dbPath, kmerFiles_c[[i]], count_matrix)
      
      # For benchmarking, create tmp file then delete
      if(runBenchmark == TRUE){
        .write_dbTable("tmp.sqlite", "tmp", count_matrix)
        dbSize <- file.info(dbPath)$size / (1024 * 1024)
        unlink("tmp.sqlite")
      } else {
        dbSize <- NA
      }

      # Bench end
      benchE <- .bench_Func(runBenchmark, klen1Remove[i], "kmer_Find", "kCnt", dbSize, "End")
      rbind(benchS, benchE)
    }
  })
  
  # Aggregate bench results
  bench_cnt <- .bench_Aggregate(runBenchmark, bench_cnt, nrow(dfTaxa), numCores, numFolds, "kmerFind")
  
  # Garbage collection
  gc()
  
  # List for kmer position files
  kmerFiles_p <- replicate(length(klen), list())
  
  # Extract kmer position data in a matrix like format for all sequences and klengths and write to db
  bench_pos <- with_progress({
    p <- progressor(along = 1:length(klen), offset = 1)
    foreach(i=1:length(klen)) %do% {

      # Message update
      p(paste0("Extracting kmer position data and writing to db", " (klen=", klen[i], ")", sep=""))
      
      # for benchmarking
      if(runBenchmark == TRUE){
        dbSize <- file.info(dbPath)$size / (1024 * 1024)
      } else {
        dbSize <- NA
      }

      # Bench st
      benchS <- .bench_Func(runBenchmark, klen[i], "kmer_Find", "kPos", dbSize, "Start")
    
      # Grab position data with c++ function for each klen, separate out taxonomy columns, separate out fold data
      pos_matrix <- kmer_posFindCpp(as.character(dfTaxa$dna), dfTaxa$taxonomy, klen[i]) %>%
        ungroup() %>%
        separate(SeqID, nameVecTax, sep = ";")
      
      # If cross-val tests are being performed, organize the fold data first into train or test
      if(numFolds >= 3){
        # Character for folds
        foldCols <- paste0("Fold", 1:numFolds, sep="")
        
        # Organize fold data for training
        pos_matrix_train <- rbindlist(foreach(j=1:numFolds) %do% {
          # Subset the data into different folds
          pos_matrix[which(pos_matrix[,j+5] == "train"),] %>%
            mutate(folds = paste0("f", j), sep="")
        }) %>% select(!all_of(foldCols))
        
        # Organize fold data for testing
        pos_matrix_test <- rbindlist(foreach(j=1:numFolds) %do% {
          # Subset the data into different folds
          pos_matrix[which(pos_matrix[,j+5] == "test"),] %>%
            mutate(folds = paste0("f", j, sep=""))
        }) %>% select(!all_of(foldCols))
        
        # Replace NA's with gaps
        pos_matrix_test[is.na(pos_matrix_test)] <- "-"
        pos_matrix_train[is.na(pos_matrix_train)] <- "-"
        
        # Write out to db for each fold, train/test etc. and store table names in a list
        foreach(j=1:numFolds) %do% {
          # write test data to db
          .write_dbTable(dbPath, 
                         paste("raw_klen_", klen[i],  "_f", j, "_test_pos", sep=""), 
                         pos_matrix_test %>% 
                           filter(folds == paste0("f", j, sep="")) %>% 
                           select(!folds))
          # write train data to db
          .write_dbTable(dbPath, 
                         paste("raw_klen_", klen[i],  "_f", j, "_train_pos", sep=""), 
                         pos_matrix_test %>% 
                           filter(folds == paste0("f", j, sep="")) %>% 
                           select(!folds))
        }
        
        # Sublists for kmer file p
        kmerFiles_p[[i]] <- replicate(numFolds, list())

        # Write out names for kPosdata tables to list
        foreach(j=1:numFolds) %do% {
          kmerFiles_p[[i]][[j]][["test"]] <- paste0("raw_klen_", klen[i],  "_f", j, "_test_pos", sep="")
          kmerFiles_p[[i]][[j]][["train"]] <- paste0("raw_klen_", klen[i],  "_f", j, "_train_pos", sep="")
        }
        names(kmerFiles_p[[i]]) <- paste0("f", 1:numFolds, sep="")
        
      # Else if not performing cross-validation tests
      } else {
        # write test to db
        .write_dbTable(dbPath, 
                       paste("raw_klen_", klen[i], "_pos", sep=""), 
                       pos_matrix)
        # store table names
        kmerFiles_p[[i]] <- paste("raw_klen_", klen[i], "_pos", sep="")
      }
      
      # For benchmarking
      if(runBenchmark == TRUE){
        dbSize <- file.info(dbPath)$size / (1024 * 1024)
      } else {
        dbSize <- NA
      }
      
      # Bench end
      benchE <- .bench_Func(runBenchmark, klen[i], "kmer_Find", "kPos", dbSize, "End")
      rbind(benchS, benchE)
    }
  })
  
  # Garbage collection
  gc()
  
  # Aggregate bench results
  bench_pos <- .bench_Aggregate(runBenchmark, bench_pos, nrow(dfTaxa), 1, numFolds, "kmer_Find")
  
  # Depending on if klen 1 is included
  if(sum(klen)!=1){
    # Name files
    names(kmerFiles_c) <- paste0("klen", klen1Remove, sep="")
    names(kmerFiles_p) <- paste0("klen", klen, sep="")
    kmerFiles <- list("CountData"=kmerFiles_c, "PosData"=kmerFiles_p)
    
    # Remove length 0 entries for filenames in count data
    kmerFiles$CountData <- kmerFiles$CountData[lapply(kmerFiles$CountData, length)>0]
  } else {
    names(kmerFiles_p) <- paste0("klen", klen, sep="")
    kmerFiles <- list("CountData"=kmerFiles_c, "PosData"=kmerFiles_p)
  }
  
  # Write bench results to db
  if(runBenchmark == TRUE){
    .write_dbTable(dbPath, "bench_kmer_Find", rbind(bench_cnt, bench_pos)) 
    kmerFiles[["BenchData"]] <- c("bench_kmer_Find")
  }
  suppressWarnings(rm(bench_pos, bench_cnt, raw_bench))
  
  # Final message
  cat(green(paste0("Kmers (klen=", paste(klen, collapse=','), ") have been extracted and written to the ", dir, 
                   ".sqlite db.\n")))
  
  # Remove any remaining vars not needed
  rm(ids, nameVecFolds, nameVecTax)
  
  # return files created for each kmer length
  return(kmerFiles)
}

# Function to find numMaxK depending on klen
calc_numMaxK <- function(klenSub, scaleCoefMaxK){
  numMaxK <- c()
  foreach(i=1:length(klenSub)) %do% {
    if(klenSub[i]==2){
      # klen2 value
      numMaxK[i] <- ceiling(length(.find_Residues(klenSub[i])) * scaleCoefMaxK[1])
    } else if(klenSub[i]==3) {
      # klen3 value
      numMaxK[i] <- ceiling(length(.find_Residues(klenSub[i])) * scaleCoefMaxK[2])
    } else if(klenSub[i]==4) {
      # klen4 value
      numMaxK[i] <- ceiling(length(.find_Residues(klenSub[i])) * scaleCoefMaxK[3])
    } else if(klenSub[i]==5) {
      # klen5 value
      numMaxK[i] <- ceiling(length(.find_Residues(klenSub[i])) * scaleCoefMaxK[4])
    } else if(klenSub[i]>5){
      # klen>5 value
      numMaxK[i] <- ceiling(length(.find_Residues(klenSub[i])) * scaleCoefMaxK[5])
    }
  }
  return(numMaxK)
}

# Main function for finding discriminative kmers - uses two internal functions .kmer_wilcoxTest & .kmer_Filter
kmer_discrimFind <- function(kmerFiles, dfTaxa, numFolds, subclassRank, 
                             klen, maxSigVal, numMinK, numMaxK, runParallel, 
                             numCores, futureGlobals, dbPath, runBenchmark){
  
  # Message update
  cat(green(paste0("Importing kmer data...\n")))
  
  # Initialize parallel processing settings
  if(runParallel == TRUE){
    options(future.globals.maxSize = futureGlobals)
    registerDoFuture()
    plan(multisession, workers=numCores)
    handlers("progress")
  } else {
    plan(sequential)
    handlers("progress")
  }
  
  # List for wilcox test results data
  wTest <- replicate(length(klen), list())
  
  # Empty list for benchmarking
  bench <- list()

  # Iterate through each kmer length extracting the results
  foreach(i=1:length(klen)) %do% {

    # Bench st
    benchS <- .bench_Func(runBenchmark, klen[i], "kmer_discrimFind", "None", 0, "Start")
    
    # Stop function if numMinK > numMaxK for any klen
    if(numMinK > numMaxK[i]){
      stop("numMinK cannot be higher than numMaxK for any klen.")
    }
    
    # Read in raw kmer data
    kData <- read_dbTable(dbPath, kmerFiles[[1]][[i]])
    
    # Set kData as a data.table
    setDT(kData)
    
    # Perform wilcox rank sum tests depending on which taxonomic ranks are being run
    # Order level
    if("Order" %in% subclassRank){
      rank <- "Order"
      
      # Find all unique orders
      uniqOrd <- sort(unique(kData$ordID))
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # Iterate through each unique order
      wTest[[i]][[rank]] <- with_progress({
        p <- progressor(along = 1:length(uniqOrd), offset = 1)
        rbindlist(foreach(j=1:length(uniqOrd)) %dofuture% { 
          # Update progress
          p(paste0("Searching kmers... klen=", klen[i], ", Rank: ", rank, sep=""))
          # Run wilcox tests
          tryCatch({ .kmer_wilcoxTest(kData[ordID == uniqOrd[j]], 
                                      kData[ordID != uniqOrd[j]], 
                                      uniqOrd[j], rank, numFolds, maxSigVal) }, error = function(e) { NULL })
        })})
      
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK[i])
    }
    
    # Family Level
    if("Family" %in% subclassRank){
      rank <- "Family"
      
      # Find all unique orders
      uniqFam <- sort(unique(kData$famID))
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # Iterate through each unique order
      wTest[[i]][[rank]] <- with_progress({
        p <- progressor(along = 1:length(uniqFam), offset = 1)
        rbindlist(foreach(j=1:length(uniqFam)) %dofuture% { 
          # Update progress
          p(paste0("Searching kmers... klen=", klen[i], ", Rank: ", rank, sep=""))
          # Run wilcox tests
          tryCatch({ .kmer_wilcoxTest(kData[famID == uniqFam[j]], 
                                      kData[famID != uniqFam[j]], 
                                      uniqFam[j], rank, numFolds, maxSigVal) }, error = function(e) { NULL })
        })})
      
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK[i])
    }
    
    # Genus Level
    if("Genus" %in% subclassRank){
      rank <- "Genus"
      
      # Find all unique orders
      uniqGen <- sort(unique(kData$genID))
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # Iterate through each unique order
      wTest[[i]][[rank]] <- with_progress({
        p <- progressor(along = 1:length(uniqGen), offset = 1)
        rbindlist(foreach(j=1:length(uniqGen)) %dofuture% { 
          # Update progress
          p(paste0("Searching kmers... klen=", klen[i], ", Rank: ", rank, sep=""))
          # Run wilcox tests
          tryCatch({ .kmer_wilcoxTest(kData[genID == uniqGen[j]], 
                                      kData[genID != uniqGen[j]], 
                                      uniqGen[j], rank, numFolds, maxSigVal) }, error = function(e) { NULL })
        })})
      
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK[i])
    }
    
    # Species Level
    if("Species" %in% subclassRank){
      rank <- "Species"
      
      # Find all unique orders
      uniqSp <- sort(unique(kData$spID))
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # Iterate through each unique order
      wTest[[i]][[rank]] <- with_progress({
        p <- progressor(along = 1:length(uniqSp), offset = 1)
        rbindlist(foreach(j=1:length(uniqSp)) %dofuture% { 
          # Update progress
          p(paste0("Searching kmers... klen=", klen[i], ", Rank: ", rank, sep=""))
          # Run wilcox tests
          tryCatch({ .kmer_wilcoxTest(kData[spID == uniqSp[j]], 
                                      kData[spID != uniqSp[j]], 
                                      uniqSp[j], rank, numFolds, maxSigVal) }, error = function(e) { NULL })
        })})
      
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK[i])
    }
    
    # Bench end
    benchE <- .bench_Func(runBenchmark, klen[i], "kmer_discrimFind", "None", 0, "End")
    bench[[i]] <- rbind(benchS, benchE)
    
    # Garbage collection
    gc()
  }
  
  # Aggregate bench results
  bench <- .bench_Aggregate(runBenchmark, bench, nrow(dfTaxa), numCores, numFolds, "kmer_discrimFind")
  
  # Rbindlist of discriminative kmer data into one dataframe
  wTest <- rbindlist(foreach(i=1:length(wTest)) %do% {
    rbindlist(wTest[[i]]) %>% mutate(klen = klen[i])
  }) 
  
  # Tally for numbers of discrim kmers per group
  if(numFolds >= 3){
    tallyK <- wTest %>%
      select(klen, group, rank, fold, kCount_Group) %>%
      distinct()
  } else {
    tallyK <- wTest %>%
      select(klen, group, rank, kCount_Group) %>%
      distinct()
  }
  
  # Groups where kmers could not be found at each rank, klen and fold if applicable
  groupDict <- rbind(
    data.frame("group" = unique(dfTaxa$ordID), "rank" = "Order"),
    data.frame("group" = unique(dfTaxa$famID), "rank" = "Family"),
    data.frame("group" = unique(dfTaxa$genID), "rank" = "Genus"),
    data.frame("group" = unique(dfTaxa$spID), "rank" = "Species")
  )
  if(numFolds >= 3){
    noK <- wTest %>% 
      select(c(group, rank, fold, klen)) %>% 
      group_by(fold, rank, klen) %>%
      group_split()
  } else {
    noK <- wTest %>% 
      select(c(group, rank, klen)) %>% 
      group_by(rank, klen) %>%
      group_split()
  }
  foreach(i=1:length(noK)) %do% {
    noK[[i]] <- groupDict %>% 
      filter(rank == noK[[i]]$rank[1]) %>% 
      filter(str_detect(group, paste0(unique(noK[[i]]$group), collapse="|"), negate=TRUE)) %>%
      mutate(klen = noK[[i]]$klen[1]) %>% 
      mutate(fold = ifelse(numFolds >= 3, noK[[i]]$fold[1], NA))
  }
  noK <- rbindlist(noK)

  # Remove columns not needed
  wTest <- wTest %>% select(!c("kCount_Group"))
  
  # kmerResults list to return from the function
  kmerResults <- list("kResults" = wTest,
                      "kTally" = tallyK,
                      "kFail" = noK)
  
  # Write bench results to db if runBenchmark is TRUE
  if(runBenchmark == TRUE){
    .write_dbTable(dbPath, "bench_kmer_discrimFind", bench) 
    kmerResults[["BenchData"]] <- c("bench_kmer_discrimFind")
  }
  suppressWarnings(rm(bench))
  
  # Remove unneeded vars
  suppressWarnings(rm(wTest, tallyK, kData, noK, groupDict))
  
  # Return kmerResults
  return(kmerResults)
}

# Internal function to perform wilcoxon ranked sum tests

# For each kmer group and taxonomic group at a specified taxonomic rank, perform a series of 
# group vs non-group comparisons. In this case group refers to kmer frequencies of target taxonomic group
# and non-group is any other kmer frequencies not belonging to that group (ex: all Bovidae kmer freqs vs all Non-Bovidae kmer freqs)
# A two-sided, non-paired wilcoxon rank-sum test is then performed for each group/nongroup and p-values are extracted from the test
# and corrected via bonferonni correction
.kmer_wilcoxTest <- function(col_g, col_ng, group, rank, numFolds, maxSigVal){
  # If cross-val tests are not being performed
  if(numFolds < 3){
    if(!is.null(ncol(col_g))){
      wilcox <- suppressWarnings(matrixTests::col_wilcoxon_twosample(
        as.matrix(col_g[,6:ncol(col_g)]),
        as.matrix(col_ng[,6:ncol(col_ng)]),
        null = 0,
        alternative = "two.sided",
        exact = NA,
        correct = TRUE)) %>% 
        rownames_to_column(var = "kmers") %>%
        adjust_pvalue(method = "bonferroni") %>%
        # Filter by maxSigVal param
        filter(pvalue.adj < maxSigVal) %>%
        mutate(group = group) %>% 
        mutate(rank = rank) %>% 
        select(group, rank, kmers, pvalue.adj)
    } else {
      wilcox <- NULL
    }
  # If cross-val tests are being performed
  } else {
    maxPos <- 5+numFolds
    positions <- c(1,1:maxPos)
    # Find the rows labelled as train from each fold and perform tests on them only
    wilcox <- foreach(i=1:numFolds) %do% {
      if(!is.null(ncol(col_g))){
        suppressWarnings(matrixTests::col_wilcoxon_twosample(
          as.matrix(col_g %>% filter(select(.,i+5) == "train") %>% select(!all_of(positions))),
          as.matrix(col_ng %>% filter(select(.,i+5) == "train") %>% select(!all_of(positions))),
          null = 0,
          alternative = "two.sided",
          exact = NA,
          correct = TRUE)) %>% 
          rownames_to_column(var = "kmers") %>%
          adjust_pvalue(method = "bonferroni") %>%
          # Filter by maxSigVal param
          filter(pvalue.adj < maxSigVal) %>%
          mutate(group = group) %>% 
          mutate(rank = rank) %>% 
          select(group, rank, kmers, pvalue.adj) %>%
          mutate(fold = paste0("f", i, sep=""))
      } else {
        wilcox <- NULL
      }
    }
    if(!is.null(wilcox)){
      # Lastly rbindlist
      wilcox <- rbindlist(wilcox)
    } else {
      wilcox <- NULL
    }
  }
  
  # remove unneeded vars
  suppressWarnings(rm(maxPos, positions))
  
  # return test results
  return(wilcox)
}

# Internal function to filter discriminative kmers by numMaxK
.kmer_Filter <- function(wTest, numFolds, numMinK, numMaxK){
  
  # Number of unique groups present
  uniqGroups <- length(unique(wTest$group))
  
  # If cross-validation tests are being performed
  if(numFolds >=3){
    # Determine kCoverage (proportion of all unique groups that possess said kmer) and filter by kCount_Group >= numMinK
    setDT(wTest)
    # Group by 'kmers' and add 'kCoverage'
    wTest[, kCoverage := .N / uniqGroups, by = c("fold", "kmers")]
    # Group by 'group' and add 'kCount_Group'
    wTest[, kCount_Group := .N, by = c("fold", "group")]
    # Filter and arrange
    wTest <- wTest[kCount_Group >= numMinK[1],][order(kCoverage, pvalue.adj)][,head(.SD, numMaxK), by = group]
    # Recalculate numbers per group and fold 
    wTest[, kCount_Group := .N, by = c("fold", "group")]
    
  # If no cross-validation tests are not being performed
  } else {
    # Determine kCoverage (proportion of all unique groups that possess said kmer) and filter by kCount_Group >= numMinK
    setDT(wTest)
    # Group by 'kmers' and add 'kCoverage'
    wTest[, kCoverage := .N / uniqGroups, by = kmers]
    # Group by 'group' and add 'kCount_Group'
    wTest[, kCount_Group := .N, by = group]
    # Filter and arrange
    wTest <- wTest[kCount_Group >= numMinK[1],][order(kCoverage, pvalue.adj)][,head(.SD, numMaxK), by = group]
    # Recalculate numbers per group and fold
    wTest[, kCount_Group := .N, by = group]
  }
  
  # remove vars not needed
  rm(uniqGroups)
  
  # Return wTest
  return(wTest)
}

##### kPHHM Model Setup and Training Functions #####

# Main function used to do all setup & training of kPHMM models
kphmm_Train <- function(kmerResults, kmerFiles, dfTaxa, klen, numFolds, numMaxiter, 
                        subclassRank, gapThreshold, deltaLL, bepK, runParallel, numCores, 
                        futureGlobals, doRNGseed, dir, dbPath, writeHMM, mixedModels, 
                        numSampleTrain, runBenchmark){

  # Message update
  cat(green(paste0("Importing kmer data...\n")))
  
  # Initialize parallel processing settings
  if(runParallel == TRUE){
    options(future.globals.maxSize = futureGlobals)
    registerDoFuture()
    plan(multisession, workers=numCores)
    handlers("progress")
  } else {
    plan(sequential)
    handlers("progress")
  }
  
  # Creating directory for model .HMM files if writeHMM is TRUE
  if(writeHMM == TRUE){
    dirM <- paste0(dir, "_ModelData")
    dir.create(dirM, showWarnings = FALSE)
  } else {
    dirM <- NA
  }
  
  # Grab kmer Results data and split by rank & group, organizing in the correct taxonomic rank order
  # from subclassRank
  kmers <- kmerResults[['kResults']] %>%
    mutate(rank = factor(rank, levels = subclassRank, 
                         labels = subclassRank)) %>%
    group_by(rank, group) %>%
    group_split()
  
  # For random sampling of taxonomic groups with sampleNumTrain param if it is not NA
  if(!is.na(numSampleTrain)){
    kmers <- sample(kmers, numSampleTrain)
  }
  
  # List for residues
  residues <- replicate(length(klen), list())
  
  # Creating a list of residues for use in the models
  if(1 %in% klen){
    # If klen 1 included in k lengths, it uses the "DNA" argument instead
    # of custom residues
    residues[[1]] <- "DNA"
    names(residues)[1] <- "1"
    if(sum(klen)>1){
      # Create residues for all klens to be used in the PHMMs
      foreach(i=2:length(klen)) %do% {
        residues[[i]] <- .find_Residues(klen[i])
        names(residues)[i] <- klen[i]
      }
      # If mixModels is TRUE, add another residues list element,
      # representing residues of all kmer lengths (except klen 1) combined together
      if(mixedModels==TRUE){
        residues[["MixedKlen"]] <- unlist(unname(residues[2:length(residues)]))
      }
    }
  } else {
    # Create residues for all klens to be used in the PHMMs
    foreach(i=1:length(klen)) %do% {
      residues[[i]] <- .find_Residues(klen[i])
    }
    names(residues) <- klen
    # If mixModels is TRUE, add another residues list element,
    # representing residues of all kmer lengths combined together
    if(mixedModels==TRUE){
      residues[["MixedKlen"]] <- unlist(unname(residues[2:length(residues)]))
    }
  }
  
  # Message update
  #cat("\n")
  cat(green(paste0("Starting model training...\n")))
  
  # Outermost list names
  lnames <- unlist(foreach(i=1:length(kmers)) %do% paste0(kmers[[i]]$rank[1], ";", kmers[[i]]$group[1], sep=""))
  
  # Subset kmer position data for each group, klen, fold etc. generate identity matrices and then train kphmm's on those matrices
  with_progress({
    p <- progressor(along = 1:length(kmers), offset = 1)
    
    # Empty list for training results (updated globally)
    kModels <<- replicate(length(kmers), list())
    
    # Benchmark results
    if(runBenchmark == TRUE){
      benchData <<- replicate(length(kmers), list())
    }
    
    # Generate models, iterating over each set of discrim kmers for each group, klen and fold
    foreach(i=1:length(kmers)) %do% {
      
      # Update progress
      p(paste("Training kPHMMs... Rank: ", kmers[[i]]$rank[1], ", Group: ", kmers[[i]]$group[1], sep=""))
      
      # Set DT for kmers
      setDT(kmers[[i]])
      
      # Group name 
      gDat <- kmers[[i]]$group[1]
      
      # Rank name
      rDat <- kmers[[i]]$rank[1]
      
      # Unique Folds
      kFolds <- sort(unique(kmers[[i]]$fold))
      
      # If performing cross-validations
      if(numFolds >= 3){
        
        # For each fold - using dorng for parallel processing of training across folds
        kModels[[i]] <- setNames(foreach(j=1:length(kFolds), .options.RNG = doRNGseed) %dorng% {
          
          # Unique kmer lengths
          klenDat <- unique(kmers[[i]][(fold == kFolds[j])]$klen)
          
          # Add klen 1 if it is included
          if(1 %in% klen){
            klenDat <- append(1, klenDat)
          }
          
          # For each klen 
          setNames(foreach(l=1:length(klenDat)) %do% {
            
            # Bench st
            benchS <- .bench_Func(runBenchmark, klenDat[l], "kphmm_Train", "None", 0, "Start")
            
            # Table name from db
            tableNameIn <- paste0("raw_klen_", klenDat[l], "_", kFolds[j], "_train_pos", sep="")
            
            # Read in position data from db and setDT
            kPos <- read_dbTable(dbPath, tableNameIn)
            setDT(kPos)
            
            # Extract out relevant group data
            kPos <- kPos[(ordID == gDat | famID == gDat | genID == gDat | spID == gDat)]
            
            # FoldName 
            fDat <- kFolds[j]
            
            # If klen does not equal 1, subset by columns that contain discrim kmers, 
            # convert to matrix format, send recordID to rowname
            if(klenDat[l] != 1){
              # Unique kmers
              kDat <- kmers[[i]][(klen == klenDat[l] & fold == kFolds[j])]$kmers
              
              # Create a pattern that matches any of the strings
              kPattern <- paste(kDat, collapse = "|")
              
              # Find columns that contain the kmer search pattern in at least one row and filter by those columns
              kCols <- names(kPos)[sapply(kPos, function(col) any(grepl(kPattern, col)))]
              kCols <- kCols[!kCols %in% c(1,2,3,4,5)]
              
              # If number of column matches is greater than 0, convert to matrix format and send recordID
              # to rownname
              if(length(kCols)!=0){
                kPos <- cbind(kPos[,1], kPos[, ..kCols])
                kPos <- as.matrix(as.data.frame(kPos) %>% column_to_rownames(var="recordID"))
              # Else if nothing is found in kCols, assign NULL to kPos
              } else {
                kPos <- NULL
              }
              
            # For klen 1 convert to matrix format, send recordID to rowname, assign NA to kDat and kPattern 
            } else {
              kPos <- kPos[, c("ordID", "famID", "genID", "spID") := NULL]
              kPos <- as.matrix(as.data.frame(kPos) %>% column_to_rownames(var="recordID"))
              kDat <- NA
              kPattern <- NA
            }
            
            # If kPos is null, assign False for model training success
            if(is.null(kPos)){
              kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                            "klen" = klenDat[l], "fold" = fDat, "model" = NULL, 
                            "success" = FALSE)
              
            # Else if kPos ncol + nrow < 2, assign False for model training success
            } else {
              if((ncol(kPos) + nrow(kPos))<2){
                kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                              "klen" = klenDat[l], "fold" = fDat, "model" = NULL, 
                              "success" = FALSE)
                
              # Else begin training
              } else {
                # Setting numcores to 1 here since the fold loop already uses parallel processing
                kphmm <- .kphmm_runTrain(kPos, kDat, kPattern, klenDat[l], gDat, rDat, fDat,
                                         unname(unlist(residues[names(residues)==klenDat[l]])),
                                         numFolds, gapThreshold, numMaxiter, deltaLL, 
                                         bepK, dirM, dbPath, writeHMM, 1)

              }
            }
            
            # Bench End, only record if success is TRUE for the model, else ignore
            if(kphmm$success == TRUE){
              modelNameOut <- paste0("model_", gDat, "_", fDat, "_", klenDat[l], sep="")
              benchE <- .bench_Func(runBenchmark, klenDat[l], "kphmm_Train", "None", 0, "End")
              bench <- rbind(benchS, benchE)
            } else {
              modelNameOut <- NA
              bench <- NULL
            }
            
            # model name, model and bench data as last line
            if(runBenchmark == TRUE){
              list("modelName" = modelNameOut,
                   "modelData" = kphmm,
                   "benchData" = bench)
            } else {
              list("modelName" = modelNameOut,
                   "modelData" = kphmm)
            }
            
            # Name with klens
          }, klenDat)
          
          # Name with folds
        }, kFolds)
        
        # If the mixed models option is true, generate mixed models and add to each individual fold
        if(mixedModels == TRUE){
          # For each fold
          kMixModels <- setNames(foreach(j=1:length(kFolds), .options.RNG = doRNGseed) %dorng% {
            
            # Start benchmark
            benchS <- .bench_Func(runBenchmark, "MixedKlen", "kphmm_Train", "None", 0, "Start")
            
            # For fold data
            fDat <- kFolds[j]
            
            # remove klen=1 if it is present
            matMix <- kModels[[i]][[j]][names(kModels[[i]][[j]])!=1]
            
            # Unique set of klengths without 1
            uniqKlenMix <- names(matMix)
            
            # Extract out mats from matMix
            matMix <- foreach(l=1:length(matMix)) %do% as.data.frame(matMix[[l]]$modelData$model$alignment)
            
            # Unique set of klengths without 1
            names(matMix) <- uniqKlenMix
            
            # Filter out NA's
            matMix <- matMix[lapply(matMix, nrow)>1]
            
            # Store klen's used
            klenDatMixed <- names(matMix)
            
            # If only one kmer length (besides klen 1), skip this step
            if(length(matMix)>1){
              # Modify the column names of each klen
              foreach(l=1:length(matMix)) %do% {
                colnames(matMix[[l]])[1:ncol(matMix[[l]])] <- paste0(uniqKlenMix[l], "_", 
                                                                     colnames(matMix[[l]][,1:ncol(matMix[[l]])]), sep="")
              }
              
              # Col bind dfs together
              matMix <- dplyr::bind_cols(matMix)
              
              # Name for matMix
              modelName <- paste0("model_", gDat, "_", fDat, "_MixedKlen", sep="")
              
              # Unique discrim kmers used by the mixed klen model
              kDat <- kmers[[i]][(fold == fDat)]$kmers
              kPattern <- paste0(kDat, collapse="|")
              
              # if ncol + nrow < 2, assign False for model training success
              if((ncol(matMix) + nrow(matMix))<2){
                kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                              "klen" = "MixedKlen", "fold" = fDat, "model" = NULL, "success" = FALSE)
                bench <- NULL
                modelName <- NA
                
              # Else begin training  
              } else {
                # Train a mixed klen model, use the numCores option for parallelization in this step
                kphmm <- .kphmm_runTrain(matMix, kDat, kPattern, "MixedKlen", gDat, rDat, fDat,
                                         unname(unlist(residues[names(residues)=="MixedKlen"])),
                                         numFolds, gapThreshold, numMaxiter, deltaLL, bepK, 
                                         dirM, dbPath, writeHMM, numCores)
                
                # If training was successful
                if(kphmm$success == TRUE){
                  # For benchmarking
                  if(runBenchmark == TRUE){
                    benchE <- .bench_Func(runBenchmark, "MixedKlen", "kphmm_Train", "None", 0, "End")
                    bench <- rbind(benchS, benchE)
                  } else {
                    bench <- NULL
                  }
                  
                # If training wasn't successful
                } else {
                  bench <- NULL
                  modelName <- NA
                }
              }
              
              # Create list with model details including bench data if runBenchmark is TRUE
              if(runBenchmark == TRUE){
                list("modelName"= modelName, "modelData" = kphmm, "benchData" = bench)
              } else {
                list("modelName"= modelName, "modelData" = kphmm)
              }
              
            # Else do nothing, NULL will be assigned
            } else {
            }
            
          # Name with folds
          }, kFolds)
          
          # Add the mixed models to the kModels list
          # kModels <- setNames(foreach(j=1:length(kModels)) %do% {
          foreach(j=1:length(kModels[[i]])) %do% {
            if(!is.null(kMixModels[[j]])){
              kModels[[i]][[j]][["MixedKlen"]] <- kMixModels[[j]]
            }
          }
          
          # Remove mixed models
          rm(kMixModels)
        }
        
        # Garbage collection
        gc()
        
        # If benchmarking is TRUE run some tests to quantify db size
        if(runBenchmark == TRUE){
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            
            # Naming var for innermost loop
            klenDat <- names(kModels[[i]][[j]])
            
            # For each klen
            setNames(foreach(l=1:length(kModels[[i]][[j]])) %do% {
              if(kModels[[i]][[j]][[l]]$modelData$success == TRUE){
                # For benchmark of db size before write and then append dbSize to benchData
                dbSizeSt <- file.info(dbPath)$size / (1024 * 1024)
                
                # write out to db
                .write_dbTable(dbPath, 
                               kModels[[i]][[j]][[l]]$modelName, 
                               kModels[[i]][[j]][[l]]$modelData$model)
                
                # For benchmark of db size after write and then append dbSize to benchData
                dbSizeEnd <- file.info(dbPath)$size / (1024 * 1024)
                kModels[[i]][[j]][[l]]$benchData[stORend=="End"]$dbSize <- dbSizeEnd - dbSizeSt
                
                # Replace the model itself with the name of the table written to db
                # (since its too large to keep in memory)
                kModels[[i]][[j]][[l]]$modelData$model <- kModels[[i]][[j]][[l]]$modelName
              }
              
              # Get rid of model name
              kModels[[i]][[j]][[l]][names(kModels[[i]][[j]][[l]])!="modelName"]
              
            # Name with klength
            }, klenDat)
            
          # Name with folds
          }, kFolds)
          
          # Aggregate bench data and send to a different list
          benchData[[i]] <- rbindlist(foreach(j=1:length(kModels[[i]])) %do% {
            rbindlist(foreach(l=1:length(kModels[[i]][[j]])) %do% {
              # Get rid of model data
              kModels[[i]][[j]][[l]][names(kModels[[i]][[j]][[l]])!="modelData"]
              kModels[[i]][[j]][[l]]$benchData
            })
          })
          
          # Name bench data with rank and group
          names(benchData)[i] <- lnames[[i]]
          
          # Remove bench data for kModels list
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            
            # Naming var for innermost loop
            klenDat <- names(kModels[[i]][[j]])
            
            # For each klen
            setNames(foreach(l=1:length(kModels[[i]][[j]])) %do% {
              
              # Get rid of model name
              kModels[[i]][[j]][[l]][names(kModels[[i]][[j]][[l]])!="benchData"]
              
            # Name with klength
            }, klenDat)
            
          # Name with folds  
          }, kFolds)
          
        # Else just write to db and remove model name
        } else {
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            
            # Naming var for innermost loop
            klenDat <- names(kModels[[i]][[j]])
            
            # For each klen
            setNames(foreach(l=1:length(kModels[[i]][[j]])) %do% {
              if(kModels[[i]][[j]][[l]]$modelData$success == TRUE){
                # write out to db
                .write_dbTable(dbPath, 
                               kModels[[i]][[j]][[l]]$modelName, 
                               kModels[[i]][[j]][[l]]$modelData$model)
                # Replace the model itself with the name of table for the model data written to db
                # (since its too large to keep in memory)
                kModels[[i]][[j]][[l]]$modelData$model <- kModels[[i]][[j]][[l]]$modelName
              }
              # Get rid of model name
              kModels[[i]][[j]][[l]][names(kModels[[i]][[j]][[l]])!="modelName"]
              
            # Name with klength  
            }, klenDat)
            
          # Name with folds              
          }, kFolds)
        }   
        
      # If not performing cross-validations (parallel processing is performed at the level of kmer length iteration instead)
      } else {
        # For each klen using dorng for parallel processing of training across individual klengths
        kModels[[i]] <- setNames(foreach(l=1:length(klenDat), .options.RNG = doRNGseed) %dorng% {
          
          # Bench st
          benchS <- .bench_Func(runBenchmark, klenDat[l], "kphmm_Train", "None", 0, "Start")
          
          # Table name from db
          tableNameIn <- paste0("raw_klen_", klenDat[l], "_pos", sep="")
          
          # Read in position data from db and setDT
          kPos <- read_dbTable(dbPath, tableNameIn)
          setDT(kPos)
          
          # Extract out relevant group data
          kPos <- kPos[(ordID == gDat | famID == gDat | genID == gDat | spID == gDat)]
          
          # If klen does not equal 1, subset by columns that contain discrim kmers, 
          # convert to matrix format, send recordID to rowname
          if(klenDat[l] != 1){
            
            # Unique kmers
            kDat <- kmers[[i]][(klen == klenDat[l])]$kmers
            
            # Create a pattern that matches any of the strings
            kPattern <- paste(kDat, collapse = "|")
            
            # Find columns that contain the kmer search pattern in at least one row and filter by those columns
            kCols <- names(kPos)[sapply(kPos, function(col) any(grepl(kPattern, col)))]
            kCols <- kCols[!kCols %in% c(1,2,3,4,5)]
            
            # If number of column matches is greater than 0, convert to matrix format and send recordID
            # to rownname
            if(length(kCols)!=0){
              kPos <- cbind(kPos[,1], kPos[, ..kCols])
              kPos <- as.matrix(as.data.frame(kPos) %>% column_to_rownames(var="recordID"))
            # Else if nothing is found in kCols, assign NULL to kPos
            } else {
              kPos <- NULL
            }
            
          # For klen 1 convert to matrix format, send recordID to rowname, assign NA to kDat and kPattern 
          } else {
            kPos <- kPos[, c("ordID", "famID", "genID", "spID") := NULL]
            kPos <- as.matrix(as.data.frame(kPos) %>% column_to_rownames(var="recordID"))
            kDat <- NA
            kPattern <- NA
          }
          
          # If kPos is null, assign False for model training success
          if(is.null(kPos)){
            kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                          "klen" = klenDat[l], "model" = NULL, "success" = FALSE)
            
          # Else if kPos ncol + nrow < 2, assign False for model training success
          } else {
            if((ncol(kPos) + nrow(kPos))<2){
              kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                            "klen" = klenDat[l], "model" = NULL, "success" = FALSE)
              
            # Else begin training
            } else {
              # Assign NA for fDat
              fDat <- NA
              # Setting numcores to 1 here since the fold loop already uses parallel processing
              kphmm <- .kphmm_runTrain(kPos, kDat, kPattern, klenDat[l], gDat, rDat, fDat,
                                       unname(unlist(residues[names(residues)==klenDat[l]])),
                                       numFolds, gapThreshold, numMaxiter, deltaLL, 
                                       bepK, dirM, dbPath, writeHMM, 1)
            }
          }
          
          # Bench End, only record if success is TRUE for the model, else ignore
          if(kphmm$success == TRUE){
            modelNameOut <- paste0("model_", gDat, "_", klenDat[l], sep="")
            benchE <- .bench_Func(runBenchmark, klenDat[l], "kphmm_Train", "None", 0, "End")
            bench <- rbind(benchS, benchE)
          } else {
            modelNameOut <- NA
            bench <- NULL
          }
          
          # model name, model and bench data as last line
          if(runBenchmark == TRUE){
            list("modelName" = modelNameOut,
                 "modelData" = kphmm,
                 "benchData" = bench)
          } else {
            list("modelName" = modelNameOut,
                 "modelData" = kphmm)
          }
          
        # Name with klens
        }, klenDat)
        
        # If the mixed models option is true, generate mixed models and add to each individual fold
        if(mixedModels == TRUE){
          # Start benchmark
          benchS <- .bench_Func(runBenchmark, "MixedKlen", "kphmm_Train", "None", 0, "Start")
          
          # For fold data
          fDat <- NA
          
          # remove klen=1 if it is present
          matMix <- kModels[[i]][names(kModels[[i]])!=1]
          
          # Unique set of klengths without 1
          uniqKlenMix <- names(matMix)
          
          # Extract out mats from matMix
          matMix <- foreach(l=1:length(matMix)) %do% as.data.frame(matMix[[l]]$modelData$model$alignment)
          
          # Unique set of klengths without 1
          names(matMix) <- uniqKlenMix
          
          # Filter out NA's
          matMix <- matMix[lapply(matMix, nrow)>1]
          
          # Store klen's used
          klenDatMixed <- names(matMix)
          
          # If only one kmer length (besides klen 1), skip this step
          if(length(matMix)>1){
            
            # Modify the column names of each klen
            foreach(l=1:length(matMix)) %do% {
              colnames(matMix[[l]])[1:ncol(matMix[[l]])] <- paste0(uniqKlenMix[l], "_", 
                                                                   colnames(matMix[[l]][,1:ncol(matMix[[l]])]), sep="")}
            # Col bind dfs together
            matMix <- dplyr::bind_cols(matMix)
            
            # Name for matMix
            modelName <- paste0("model_", gDat, "_MixedKlen", sep="")
            
            # Unique discrim kmers used by the mixed klen model
            kDat <- kmers[[i]]$kmers
            kPattern <- paste0(kDat, collapse="|")
            
            # if ncol + nrow < 2, assign False for model training success
            if((ncol(matMix) + nrow(matMix))<2){
              kphmm <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                            "klen" = "MixedKlen", "model" = NULL, "success" = FALSE)
              bench <- NULL
              modelName <- NA
              
            # Else begin training  
            } else {
              # Train a mixed klen model, use the numCores option for parallelization in this step
              kphmm <- .kphmm_runTrain(matMix, kDat, kPattern, "MixedKlen", gDat, rDat, fDat,
                                       unname(unlist(residues[names(residues)=="MixedKlen"])),
                                       numFolds, gapThreshold, numMaxiter, deltaLL, bepK, 
                                       dirM, dbPath, writeHMM, numCores)
              # If training was successful
              if(kphmm$success == TRUE){
                
                # For benchmarking
                if(runBenchmark == TRUE){
                  benchE <- .bench_Func(runBenchmark, "MixedKlen", "kphmm_Train", "None", 0, "End")
                  bench <- rbind(benchS, benchE)
                } else {
                  bench <- NULL
                }
                
              # If training wasn't successful
              } else {
                bench <- NULL
                modelName <- NA
              }
            }
            
            # Create list with model details including bench data if runBenchmark is TRUE
            if(runBenchmark == TRUE){
              kMixModels[[i]] <- list("modelName"= modelName, "modelData" = kphmm, "benchData" = bench)
            } else {
              kMixModels[[i]] <- list("modelName"= modelName, "modelData" = kphmm)
            }
            
          # Else assign NULL
          } else {
            kMixModels[[i]] <- NULL
          }
          
          # Add the mixed models to the kModels list
          # kModels <- setNames(foreach(j=1:length(kModels)) %do% {
          if(!is.null(kMixModels)){
            kModels[[i]][["MixedKlen"]] <- kMixModels[[i]]
          }
          
          # Remove mixed models
          rm(kMixModels)
        }
        
        # If benchmarking is TRUE run some tests to quantify db size
        if(runBenchmark == TRUE){
          
          # Naming var for innermost loop
          klenDat <- names(kModels[[i]])
          
          # For each klen
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            if(kModels[[i]][[j]]$modelData$success == TRUE){
              
              # For benchmark of db size before write and then append dbSize to benchData
              dbSizeSt <- file.info(dbPath)$size / (1024 * 1024)
              
              # write out to db
              .write_dbTable(dbPath, 
                             kModels[[i]][[j]]$modelName, 
                             kModels[[i]][[j]]$modelData$model)
              
              # For benchmark of db size after write and then append dbSize to benchData
              dbSizeEnd <- file.info(dbPath)$size / (1024 * 1024)
              kModels[[i]][[j]]$benchData[stORend=="End"]$dbSize <- dbSizeEnd - dbSizeSt
              
              # Replace the model itself with the name of table name for the model data written to db
              # (since its too large to keep in memory)
              kModels[[i]][[j]]$modelData$model <- kModels[[i]][[j]]$modelName
            }
            
            # Get rid of model name
            kModels[[i]][[j]][names(kModels[[i]][[j]])!="modelName"]
          
          # Name with klens
          }, klenDat)
          
          # Aggregate bench data and send to a different list
          benchData[[i]] <- rbindlist(foreach(j=1:length(kModels[[i]])) %do% {
            # Get rid of model data
            kModels[[i]][[j]][names(kModels[[i]][[j]])!="modelData"]
            kModels[[i]][[j]]$benchData
          })
          
          # Name bench data with rank and group
          names(benchData)[i] <- lnames[[i]]
          
          # Remove bench data for kModels list
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            
            # Get rid of model name
            kModels[[i]][[j]][names(kModels[[i]][[j]])!="benchData"]
            
          # Name with klens
          }, klenDat)
          
        # Else just write to db and remove model name
        } else {
          # Naming var for innermost loop
          klenDat <- names(kModels[[i]])
          
          # For each klen
          kModels[[i]] <- setNames(foreach(j=1:length(kModels[[i]])) %do% {
            if(kModels[[i]][[j]]$modelData$success == TRUE){
              
              # write out to db
              .write_dbTable(dbPath, 
                             kModels[[i]][[j]]$modelName, 
                             kModels[[i]][[j]]$modelData$model)
              
              # Replace the model itself with the name of table name for the model data written to db
              # (since its too large to keep in memory)
              kModels[[i]][[j]]$modelData$model <- kModels[[i]][[j]]$modelName
            }
            
            # Get rid of model name
            kModels[[i]][[j]][names(kModels[[i]][[j]])!="modelName"]
            
          # Name with klens
          }, klenDat)
        }         
      }
      # Update kModels and benchData as global in case of crashes
      .GlobalEnv$kModels[[i]] <- kModels[[i]]
      if(runBenchmark == TRUE){
        .GlobalEnv$benchData[[i]] <- benchData[[i]]
      }
    }
    # Garbage collection
    gc()
  })
  
  # Name all completed models with rank and group names
  names(kModels) <<- lnames
  
  # Message update
  cat("\n")
  cat(green(paste0("Model training complete, organizing model data...\n")))
  
  # Make a data frame for all models where discrim kmer data was either unavailable or 
  # where models failed during training
  ids <- c("rank", "group")
  modelFail <- rbindlist(foreach(i=1:length(kModels)) %do% {
    rbindlist(foreach(j=1:length(kModels[[i]])) %do% {
      rbindlist(foreach(l=1:length(kModels[[i]][[j]])) %do% {
        if(kModels[[i]][[j]][[l]]$modelData$success == FALSE){
          data.frame("rank_group" = names(kModels)[i],
                     "fold"=names(kModels[[i]])[j], 
                     "klen"=names(kModels[[i]][[j]])[l])
        } else {
          data.frame("rank_group" = NA,
                     "fold" = NA, 
                     "klen" = NA)
        }
      })
    })
  }) %>% filter(!is.na(rank_group))
  
  # Add models into one list
  modelData <- list("Models"=kModels)
  
  # If model fail had at least one row, add to modelData list
  if(nrow(modelFail)>0){
    modelData[["ModelFail"]] <- modelFail %>% separate(rank_group, ids, sep=";")
  } else {
    modelData[["ModelFail"]] <- "No models have failed during training."
  }

  # Write bench results to db if runBenchmark is TRUE
  if(runBenchmark == TRUE){
    # Aggregate benchmark data
    benchData <- .bench_Aggregate(runBenchmark, benchData, nrow(dfTaxa), numCores, numFolds, "kphmm_Train")
    # Write to db
    .write_dbTable(dbPath, "bench_kphmm_Train", benchData) 
    # Assign fully aggregated benchmark data table name to kModels
    modelData[["BenchData"]] <- c("bench_kphmm_Train")
    # rm bench data
    rm(benchData, pos = ".GlobalEnv")
  }

  # Message update
  cat("\n")
  cat(green(paste0("All models have been written to the ", dbPath, " db.\n")))
  
  # If writeHMM is TRUE
  if(writeHMM == TRUE){
    cat(green(paste0("All model files in .HMM format have been written to the ", dirM, " directory.\n")))
  }
  
  # Remove unneeded vars
  suppressWarnings(rm(kmers, lnames, kModels, modelFail, residues, pos = ".GlobalEnv"))
  gc()
  
  # Return the modelList
  return(modelData)
}

# Function to create vectors of residues for for each klen (all possible kmers for a given klen using gtools::permutations package)
.find_Residues <- function(klen){
  # Empty list for residues
  residues <- list()
  bases <- c("A","T","C","G")
  # Use the gtools package to compute all possible permutations of kmers at each a specified klen
  possibleBases <- data.frame(gtools::permutations(n = length(bases), v = bases, r = klen, repeats.allowed = T))
  colsBases <- colnames(possibleBases)
  possibleBases$seq <- apply( possibleBases[ , colsBases ] , 1 , paste , collapse = "" )
  residues <- as.character(possibleBases$seq)
  
  # remove unecessary vars
  rm(possibleBases, colsBases, bases)
  
  # return residues
  return(residues)
}

# Internal function to actually setup and train the kPHMM itself
.kphmm_runTrain <- function(kMatTr, kDat, kPattern, klenModel, gDat, rDat, fDat, residues, numFolds, 
                            gapThreshold, numMaxiter, deltaLL, bepK, dirM, dbPath, writeHMM, numCores){
  
  # kMatSetup df used for training
  kMatSetup <- kMatTr
  
  # If klenModel does not equal 1
  if(klenModel != 1){
    
    # Assigning weights to each recordID of each matrix based on their total number of discrim kmers
    # ex: a sequence with the all possible positions would be assigned the highest possible weight proportional to the total number of sequences
    # note: average weight across all sequences must equal 1 and sum of all weights must equal total number of sequences
    
    # Dataframe used to determine weights
    kCount <- kCount <- as.data.table(kMatTr, keep.rownames = "recordID")
    kCountW_Setup <- kCount[(recordID %in% rownames(kMatSetup))]
    
    # Perform rowwise count of k-mers that match the pattern, retaining recordID
    kCountW_Setup <- kCountW_Setup[, .(
      recordID,
      count = apply(.SD, 1, function(row) sum(str_detect(row, kPattern)))
    ), .SDcols = -"recordID"]
    
    # Calculate proportional counts
    kCountW_Setup <- kCountW_Setup[, .(countProp = count / max(kCountW_Setup$count, na.rm = TRUE)), by = c("recordID")]
    
    # Sum of counts for refinement
    sumCount_Setup <- sum(kCountW_Setup$count)
    
    # Sequence weights for model setup
    weights_set <- (kCountW_Setup$count / sumCount_Setup) * nrow(kCountW_Setup)
    
    # For kmer based models (klen>1) perform a model fine tuning step for kmer based models, 
    # weighting emission probabilities of discrim kmers more heavily over emission probabilities of non-discrim kmers
    
    # Assign weights to the residues representing the discrim kmers
    resD <- kDat
    resWeightD <- replicate(length(resD), bepK)
    
    # Subset out discrim kmers from the rest of the residues which are not discriminative
    resND <- residues[!residues %in% resD]
    
    # Assign weights to the other non discriminative residues 
    resWeightND <- replicate(length(resND), 1)
    
    # Recombine weights together and normalize weights into a set of probabilities
    # to get the background emission probabilities
    qe <- unlist(append(resWeightD, resWeightND))
    qe <- qe / sum(qe)
    
    # All residues recombined (discrim and non-discrim)
    residues <- append(resD, resND)
    
    # Name emission probabilities with residues
    names(qe) <- residues
    residues <- unname(residues)
    
    # Remove unneeded vars
    suppressWarnings(rm(kCountM_Setup, kCountW_Setup, 
                        sumCount_Setup, kCount, 
                        numRecordsSetup, numRecordsSetup))
  }
  
  # NULL construct.PHMM
  construct.PHMM <- NULL
  
  # Model setup
  # If running kmer length 1
  if(klenModel == 1){
    # "DNA" option used for residues for klen=1 and seqweights are always 1 for klen == 1 and qe is set to NULL - background emission probabilities are derived 
    # from the model itself
    try(construct.PHMM <- aphid::derivePHMM(kMatSetup, seqweights = NULL, residues = residues, qe = NULL, inserts = "threshold",
                                            threshold = gapThreshold, name = gDat, alignment = TRUE, progressive = TRUE, refine = "Baumwelch", 
                                            pseudocounts = "background", maxiter = numMaxiter, deltaLL = deltaLL, cores = numCores, quiet = FALSE), 
        silent = TRUE)
    
  # If running any kmer length besides length 1
  } else {
    # Custom residues inserted for kmers, custom seqweights and custom qe's are used and determined above
    try(construct.PHMM <- aphid::derivePHMM(kMatSetup, seqweights = weights_set, residues = residues, qe = qe, inserts = "threshold", 
                                            threshold = gapThreshold, name = gDat, alignment = TRUE, progressive = TRUE, refine = "Baumwelch", 
                                            pseudocounts = "background", maxiter = numMaxiter, deltaLL = deltaLL, cores = numCores, quiet = TRUE), 
        silent = TRUE)
  }
  
  # If performing cross-validation tests
  if(numFolds >= 3){
    # If running any kmer length besides length 1
    if(klenModel != 1){
      construct.PHMM_List <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = kDat,
                                  "klen" = klenModel, "fold" = fDat, "model" = NULL, "success" = NULL)
    } else {
      construct.PHMM_List <- list("Group" = gDat, "Rank" = rDat, "DiscrimKmers" = "N/A",
                                  "klen" = klenModel, "fold" = fDat, "model" = NULL, "success" = NULL)
    }
  # Else if no cross-validation
  } else {
    # If running any kmer length besides length 1
    if(klenModel != 1){
      construct.PHMM_List <- list("Group" = gDat, "Rank"= rDat, "DiscrimKmers" = kDat,
                                  "klen" = klenModel, "model" = NULL, "success" = NULL)
    } else {
      construct.PHMM_List <- list("Group" = gDat, "Rank"= rDat, "DiscrimKmers" = "N/A",
                                  "klen" = klenModel, "model" = NULL, "success" = NULL)
    }
  }
  
  # If construct.PHMM is still NULL, assign FALSE to success, else assign TRUE 
  if(is.null(construct.PHMM)){
    construct.PHMM_List$success <- "FALSE"
  } else {
    construct.PHMM_List$success <- "TRUE"
    
    # Assign model to construct.PHMM_List
    construct.PHMM_List$model <- construct.PHMM
    
    # if writeHMM param is true, write out to hmm file as well
    if(writeHMM == TRUE){
      # Write out to HMM format
      if(numFolds >= 3){
        modelFileHMM <-  paste(dirM, "/model_", 
                               gDat, "_", 
                               fDat, "_",
                               klenModel, ".hmm", sep="")
      } else {
        modelFileHMM <-  paste(dirM, "/model_", 
                               gDat, "_", 
                               klenModel, ".hmm", sep="")
      }
      # write out each kphmm
      writePHMM(construct.PHMM, file = modelFileHMM)
      construct.PHMM_List$modelHMM <- modelFileHMM
    }
  }

  # Remove unneeded vars
  suppressWarnings(rm(kMatTr, kMatSetup, qe, resD, resND, resWeightD, resWeightND, 
                      fDat, gDat, rDat, kDat, construct.PHMM, modelFileHMM, kPattern))
  gc()
  
  # Return construct.PHMM
  return(construct.PHMM_List)
}

##### kPHMM Classification Functions #####

# Main function to perform classifications with cross-validation using the previously trained kphmm models
# Classification using kphmm's with cross-validation
kphmm_ClassifyCV <- function(kmerData, modelData, numFolds, subclassRank, klen, classRank, dfTaxa, 
                             runParallel, numCores, futureGlobals, numSampleTest, doRNGseed, dbPath, 
                             classStrategy, weights, runBenchmark, majorityVoting, mvWt){
  
  # Stop the function if numFolds is less than 3 (since that would mean cv tests were not being performed)
  if(numFolds < 3){
    stop("Ensure that the numFolds numeric is greater or equal to 3 to use this function. To run classifications without performing cross-validation tests,
          please use kphmm_ClassifyI instead.")
  }
  
  # Message update
  cat(green(paste0("Preparing models and test sequences for classifications...\n")))
  
  # If runParallel is set to TRUE, register doFuture plan and set # of workers - set by numCores param, 
  # lastly set a progress bar handler
  if(runParallel == TRUE){
    options(future.globals.maxSize = futureGlobals)
    registerDoFuture()
    plan(multisession, workers = numCores)
    handlers("progress")
  } else {
    plan(sequential)
    handlers("progress")
  }
  
  # Create a dictionary with recordID and taxonomicIDs that can be referenced to determine correct IDs
  # and find taxonomy information
  recordDict <- dfTaxa %>% 
    select(recordID, ordID, famID, genID, spID) %>%
    as.data.table()
  
  # List for cData (classification data)
  cData <<- list()
  
  # Read in kPosData for each klen
  kPosData <- list()
  kPosData <- rbindlist(foreach(i=1:length(kmerData[[2]])) %do% {
    rbindlist(foreach(j=1:length(kmerData[[2]][[i]])) %do% {
      # Read in pos data from db
      kPosData[[i]] <- read_dbTable(dbPath, kmerData[[2]][[i]][[j]][["test"]])
      # Get rid of taxonomic IDs since test sequences are being used, no peaking at the taxonomy allowed!
      kPosData[[i]] <- kPosData[[i]] %>% 
        select(!c("ordID","famID","genID","spID")) %>% 
        mutate(fold = names(kmerData[[2]][[i]])[j])
    }, fill=TRUE) %>% mutate(klen = gsub("klen", "", names(kmerData[[2]])[i]))
  }, fill=TRUE) %>% as.data.table()
  kPosData[is.na(kPosData)] <- "-"  
  
  # Garbage collection
  gc()
  
  # List of unique records by fold
  if(is.na(numSampleTest)){
    unique(kPosData$recordID)
  } else {
    # running a sample of numSample recordIDs
    uniqRec <- sample(unique(kPosData$recordID), numSampleTest)
  }
  
  # randomly sort recordID to get random recordIDs for testing
  uniqRec <- sample(uniqRec, length(uniqRec))
  
  # Create naming vectors for the taxonomy column of dfTaxa
  ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
  nameVecFolds <- paste0("Fold", 1:numFolds, sep="")
  nameVecTax <- c("recordID", ids, nameVecFolds)
  
  # Create new set of fold columns for dfTaxa
  dfTaxa <- dfTaxa %>%
    ungroup() %>%
    select(recordID, taxonomy, dna) %>%
    separate(taxonomy, nameVecTax, ";")
  
  # Grab the necessary fold data first before proceeding
  foldCols <- paste0("Fold", 1:numFolds, sep="")
  
  # Grab all fold data for dfTaxa 
  dfTaxa <- rbindlist(foreach(i=1:numFolds) %do% {
      # Subset the data into different folds
      dfTaxa[which(dfTaxa[,i+5] == "train"),] %>%
        mutate(folds = paste0("f", i))
  }) %>% select(!all_of(foldCols)) %>% as.data.table()

  # Remove folds column
  rm(foldCols)

  # Seqs for blast db (blaster package) used for querying of recordID to check start and stop site difference, 
  # Sampling up to a max of 50/25/10/5 sequences for each group and fold depending on rank
  groupDict <- rbind(
    # Order
    dfTaxa[, .(group = ordID, recordID, dna, rank = "Order", folds, famID)]
    [order(folds, famID), .SD[sample(.N, min(100, .N))], by = .(group, folds, famID)][, famID := NULL],
    # Family
    dfTaxa[, .(group = famID, recordID, dna, rank = "Family", folds, genID)]
    [order(folds, genID), .SD[sample(.N, min(50, .N))], by = .(group, folds, genID)][, genID := NULL],
    # Genus
    dfTaxa[, .(group = genID, recordID, dna, rank = "Genus", folds, spID)]
    [order(folds, spID), .SD[sample(.N, min(25, .N))], by = .(group, folds, spID)][, spID := NULL],
    # Species
    dfTaxa[, .(group = spID, recordID, dna, rank = "Species", folds)]
    [order(folds), .SD[sample(.N, min(10, .N))], by = .(group, folds)]
  ) 
  
  # Message update
  cat(green(paste0("Starting classifications...\n")))
  
  # Classifications for each rank (all that apply in the subclassRank var)
  # With a progress bar, iterate through each recordID and perform classification on models through each taxonomic rank
  invisible(capture.output(cData <<- with_progress({
    p <- progressor(along = 1:length(uniqRec), offset = 1)
    foreach(i = seq_along(uniqRec), .options.RNG = doRNGseed) %dorng% {
      p(paste0("Performing classifications for RecordID: ", uniqRec[i]))
      # Run the internal function .kphmm_classifySubclassCV
      tryCatch({ .kphmm_classifySubclassCV(
                    kPosData[(recordID == uniqRec[i])],
                    modelData$Models,
                    uniqRec[i],
                    recordDict,
                    subclassRank,
                    klen,
                    dfTaxa,
                    dbPath, 
                    classStrategy, 
                    weights, 
                    mixedModels,
                    runBenchmark, 
                    groupDict,
                    majorityVoting, 
                    mvWt) }, error = function(e) { NULL })
    }
  })))
  
  # Remove unnecessary vars
  rm(kPosData)
  gc()
  
  # If runbenchmark is TRUE, separate classifications from benchmark data
  if(runBenchmark == TRUE){
    # Use the bench aggregate function to aggregate benchmark data according to klen
    benchData <- rbindlist(foreach(i=1:length(cData)) %do% cData[[i]][["benchData"]])
    benchData <<- .bench_Aggregate(runBenchmark, benchData, length(cData), numCores, numFolds, "kphmm_ClassifyCV")
    
    # Write out benchData to db
    benchTableName <- "bench_kphmm_ClassifyCV"
    .write_dbTable(dbPath, benchTableName, benchData)
    
    # Rbindlist the classification results into one large dataframe and then split by rank
    cData <<- rbindlist(foreach(i=1:length(cData)) %do% cData[[i]][["cData"]]) %>%
      mutate(rank = factor(rank, levels = c("Order", "Family", "Genus", "Species"))) %>%
      group_by(rank) %>%
      group_split()
    
  # Else just rbindlist the classification results into one large dataframe and then split by rank
  } else {
    cData <<- rbindlist(foreach(i=1:length(cData)) %do% cData[[i]][["cData"]]) %>%
      mutate(rank = factor(rank, levels = c("Order", "Family", "Genus", "Species"))) %>%
      group_by(rank) %>%
      group_split()
  }

  # Subset out instances where classifications could not be performed 
  # ex: usually because no discrim kmers exist at a given klen for a certain taxonomic group
  # and thus a model could not be trained and used for classifications
  classifyFail <- list()
  foreach(i=1:length(cData)) %do% {
    classifyFail[[i]] <- cData[[i]] %>% filter(group == "No data")
    cData[[i]] <- cData[[i]] %>% filter(group != "No data")
  }
  
  # setkeys on elements of cData list (recordID) before merging to recordDict
  foreach(i=1:length(cData)) %do% {
    setDT(cData[[i]])
    setkey(cData[[i]], "recordID")
  }
  
  # Set recordUD key for recordDict df
  setkey(recordDict, "recordID")
  
  # Creating the final data structure to be returned by this function containing all results data
  modelResults <- list("CorrectIDs" = list(),
                       "CorrectID%" = list(),
                       "AccMetrics" = list("Data"=list(),"ROC-AUC"=list()),
                       "ClassifyFail" = list())
  
  # Character for all klens and combined klens
  if(mixedModels == TRUE){
    klenDat <- append(klen, "MixedKlen")
    if(majorityVoting == TRUE){
      klenDat <- append(klenDat, "MajorityVoting")
    }
  } else {
    klenDat <- klen
  }

  # Iterate through each rfunc# Iterate through each rank, determining CorrectIDs, CorrectID%, AccMetrics
  foreach(i=1:length(cData)) %do% {  
    # Merge together with recordDict to check which groups match
    modelResults[["CorrectIDs"]][[i]] <- merge.data.table(cData[[i]], recordDict)
    
    # Create an assigned taxa column based on which rank is being used (ex: if order rank, then use the OrdID column as the assigned taxa column)
    colnames(modelResults[["CorrectIDs"]][[i]])[which(colnames(modelResults[["CorrectIDs"]][[i]]) == paste0(tolower(ifelse(cData[[i]]$rank[1] != "Species", 
                                                                                                                           substr(cData[[i]]$rank[1], 1, 3), 
                                                                                                                           substr(cData[[i]]$rank[1], 1, 2))), "ID", sep=""))] <- "actualTaxa"
    
    # Assign 1 or 0 depending on whether assignedTaxa matches actual taxa from recordDict 
    modelResults[["CorrectIDs"]][[i]] <- modelResults[["CorrectIDs"]][[i]] %>%
      rowwise() %>%
      mutate(correctID = ifelse(group == actualTaxa, 1, 0)) %>%
      select(c("recordID","z_score","group","actualTaxa","fold","correctID","klen"))
    
    # Group by group and fold and determine correct ID%'s
    modelResults[["CorrectID%"]][[i]] <- modelResults[["CorrectIDs"]][[i]] %>%
      group_by(actualTaxa, fold, klen) %>%
      dplyr::summarize(`CorrectId%` = (sum(correctID == 1) / n()) * 100, n = n())
    
    # Sublists for klen
    modelResults[["AccMetrics"]][["Data"]][[i]] <- list()
    modelResults[["AccMetrics"]][["ROC-AUC"]][[i]] <- list()
    
    # For each klen, generate ROC curve data and AUC data
    foreach(j=1:length(klenDat)) %do% {
      # Filter by klen
      data <- modelResults[["CorrectIDs"]][[i]] %>% filter(klen == klenDat[j])
      
      # If greater than 0 rows 
      if(nrow(data)>0){
        
        # If both negs and pos's present
        if(sum(data$correctID) != nrow(data) & sum(data$correctID) != 0){
          
          # Then lastly calculate TP, TN, FP, FN data for generating ROC curves with the cutpointr package
          modelResults[["AccMetrics"]][["Data"]][[i]][[j]] <- cutpointr::roc(data = data, x = z_score, class = correctID,
                                                                             pos_class = 1, neg_class = 0, direction = ">=") %>%
            # Also calculate cohens_kappa and F1 score as well
            add_metric(list(cohens_kappa, F1_score))
          
          # Determine AUC values for each set of ROC curve results using cutpointr package
          modelResults[["AccMetrics"]][["ROC-AUC"]][[i]][[j]] <- round(cutpointr::auc(modelResults[["AccMetrics"]][["Data"]][[i]][[j]]), 3)
          
          # Also calculate the MCC scores using the TP,TN,FP,FN values generated previously
          modelResults[["AccMetrics"]][["Data"]][[i]][[j]] <- modelResults[["AccMetrics"]][["Data"]][[i]][[j]] %>% 
            rowwise() %>%
            mutate(MCC_score = ((tn*tp) - (fp*fn)) / sqrt((tn+fn) * (fp+tp) * (tn+fp) * (fn+tp))) %>%
            ungroup()
          
        # Else assign NULL  
        } else {
          modelResults[["AccMetrics"]][["Data"]][[i]][[j]] <- "No usable data"
          modelResults[["AccMetrics"]][["ROC-AUC"]][[i]][[j]] <- "No usable data"
        }
        
      # Else assign NULL
      } else {
        modelResults[["AccMetrics"]][["Data"]][[i]][[j]] <- "No usable data"
        modelResults[["AccMetrics"]][["ROC-AUC"]][[i]][[j]] <- "No usable data"
      }
    }
    
    # Name each sublist with klengths
    names(modelResults[["AccMetrics"]][["Data"]][[i]]) <- klenDat
    names(modelResults[["AccMetrics"]][["ROC-AUC"]][[i]]) <- klenDat
    
    # List element for failed classifications
    modelResults[["ClassifyFail"]][[i]] <- classifyFail[[i]]
  }
  
  # Name modelResults with taxonomic rank
  names(modelResults[["CorrectIDs"]]) <- subclassRank
  names(modelResults[["CorrectID%"]]) <- subclassRank
  names(modelResults[["ClassifyFail"]]) <- subclassRank
  names(modelResults[["AccMetrics"]][["Data"]]) <- subclassRank
  names(modelResults[["AccMetrics"]][["ROC-AUC"]]) <- subclassRank
  
  # Filter out 0 nrow dfs in Fault/Error
  modelResults[["ClassifyFail"]] <- modelResults[["ClassifyFail"]][lapply(modelResults[["ClassifyFail"]],nrow)>0]
  
  # Add bench results db table to cData list if runBenchmark is TRUE
  if(runBenchmark == TRUE){
    modelResults[["BenchData"]] <- benchTableName
  }
  
  # Message update
  cat(green(paste0("Sequence classifications are now complete.\n")))
  
  # Remove unnecessary vars
  suppressWarnings(rm(recordDict, cData, benchTableName))
  gc()
  
  # Return model results
  return(modelResults)
}

#kPosDataSub <- kPosData[(recordID == uniqRec[1])]
#models <- modelData$Models
#uniqRec <- uniqRec[1]

# Internal function to run the .kphmm_classifyForwardCV function below at each taxonomic rank depending on which taxonomic ranks are being run
.kphmm_classifySubclassCV <- function(kPosDataSub, models, uniqRec, recordDict, subclassRank, 
                                      klen, dfTaxa, dbPath, classStrategy, weights, mixedModels, 
                                      runBenchmark, groupDict, majorityVoting, mvWt) {
  # Empty list for results
  cData <- list()
  
  # Find all applicable folds to test
  uniqFolds <- unique(kPosDataSub$fold)
  
  # Find all applicable klens to test, include mixed klens if TRUE
  if(mixedModels == TRUE){
    klenDat <- append(klen, "MixedKlen")
  } else {
    klenDat <- klen
  }

  # Hierarchical Strategy
  if(classStrategy == "Hierarchical"){
    # Iterate through all ranks, subsetting family by order chosen, genus by family chosen, species by genus chosen
    for (rank in subclassRank) {
      cData[[rank]] <- rbindlist(lapply(uniqFolds, function(foldD) {
        rbindlist(lapply(klenDat, function(klenD) {
          kTest <- if(klenD == "MixedKlen"){
            kTest <- kPosDataSub[(fold == foldD & klen != "1")]
          } else {
            kTest <- kPosDataSub[(fold == foldD & klen == klenD)]
          }
          modelNames <- if(rank == "Order"){
            names(models)[str_detect(names(models), rank)]
          } else {
            .kphmm_subsetbyPrevRank(cData, models, recordDict, rank, uniqRec, foldD, klenD)
          }
          if(modelNames[1] == "No data"){
            .kphmm_classifyForwardCV(kTest, NULL, uniqRec, rank, 
                                     foldD, subclassRank, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          } else {
            .kphmm_classifyForwardCV(kTest, models[modelNames], uniqRec, rank, 
                                     foldD, subclassRank, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          }
        }))
      }))
    }
    
    # If run benchmark is TRUE separate out benchmark data from classification data
    if(runBenchmark == TRUE){
      # Combine the specified columns from cData into one data.table
      bench <- rbind(cData[[1]][,c(1,4,5,7:13)],
                     cData[[2]][,c(1,4,5,7:13)],
                     cData[[3]][,c(1,4,5,7:13)],
                     cData[[4]][,c(1,4,5,7:13)],
                     cData[[1]][,c(1,4,5,14:20)],
                     cData[[2]][,c(1,4,5,14:20)],
                     cData[[3]][,c(1,4,5,14:20)],
                     cData[[4]][,c(1,4,5,14:20)])[, .SD[!duplicated(.SD)]]
      
      # Add the new column rec_fold_rank and remove recordID, fold, and rank
      bench[, rec_fold_rank := paste0(recordID, "_", fold, "_", rank)]
      bench <- bench[, !c("recordID", "fold", "rank"), with = FALSE]
      
      # Susbet cData columns to remove benchmark data
      cData[["Order"]] <- cData[["Order"]][,1:6]
      cData[["Family"]] <- cData[["Family"]][,1:6]
      cData[["Genus"]] <- cData[["Genus"]][,1:6]
      cData[["Species"]] <- cData[["Species"]][,1:6]
      
      # Append bench data to cData
      cData <- list("cData"=cData, "benchData"=bench)
      rm(bench)
    } else {
      cData <- list("cData"=cData)
    }
  }
  
  # Multi-Tiered Strategy
  if(classStrategy == "Multi-Tiered"){
    subclassHigher <- c("Order", "Family")
    subclassLower <- c("Genus", "Species")
    
    # Iterate through the higher level ranks of Order and Family, running all models at these ranks
    for (rank in subclassHigher) {
      cData[[rank]] <- rbindlist(lapply(uniqFolds, function(foldD) {
        rbindlist(lapply(klenDat, function(klenD) {
          kTest <- if(klenD == "MixedKlen"){
            kTest <- kPosDataSub[(fold == foldD & klen != "1")]
          } else {
            kTest <- kPosDataSub[(fold == foldD & klen == klenD)]
          }
          modelNames <- names(models)[str_detect(names(models), rank)]
          if(modelNames[1] == "No data"){
            .kphmm_classifyForwardCV(kTest, NULL, uniqRec, rank, 
                                     foldD, subclassHigher, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          } else {
            .kphmm_classifyForwardCV(kTest, models[modelNames], uniqRec, rank, 
                                     foldD, subclassHigher, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          }
        }))
      }))
    }
    
    # If run benchmark is TRUE separate out benchmark data from classification data
    if(runBenchmark == TRUE){
      # Combine the specified columns from cData into one data.table
      benchOF <- rbind(cData[[1]][,c(1,4,5,7:13)],
                       cData[[2]][,c(1,4,5,7:13)],
                       cData[[1]][,c(1,4,5,14:20)],
                       cData[[2]][,c(1,4,5,14:20)])[, .SD[!duplicated(.SD)]]
      
      # Add the new column rec_fold_rank and remove recordID, fold, and rank
      benchOF[, rec_fold_rank := paste0(recordID, "_", fold, "_", rank)]
      benchOF <- benchOF[, !c("recordID", "fold", "rank"), with = FALSE]
      
      # Susbet cData columns to remove benchmark data
      cData[["Order"]] <- cData[["Order"]][,1:6]
      cData[["Family"]] <- cData[["Family"]][,1:6]
    }
    
    # Use weighted averaging of scores to determine best possible Order and Family taxonomy
    cData <- .kphmm_weightAvgMulti(recordDict, cData[["Order"]], cData[["Family"]], weights)
    
    # Iterate through lower level ranks hierarchically (same as the hierarchical strategy)
    for (rank in subclassLower) {
      cData[[rank]] <- rbindlist(lapply(uniqFolds, function(foldD) {
        rbindlist(lapply(klenDat, function(klenD) {
          kTest <- if(klenD == "MixedKlen"){
            kTest <- kPosDataSub[(fold == foldD & klen != "1")]
          } else {
            kTest <- kPosDataSub[(fold == foldD & klen == klenD)]
          }
          modelNames <- .kphmm_subsetbyPrevRank(cData, models, recordDict, rank, uniqRec, foldD, klenD)
          if(modelNames[1] == "No data"){
            .kphmm_classifyForwardCV(kTest, NULL, uniqRec, rank, 
                                     foldD, subclassLower, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          } else {
            .kphmm_classifyForwardCV(kTest, models[modelNames], uniqRec, rank, 
                                     foldD, subclassLower, klenD, dfTaxa, dbPath,
                                     classStrategy, runBenchmark, groupDict)
          }
        }))
      }))
    }
    
    # If run benchmark is TRUE separate out benchmark data from classification data
    if(runBenchmark == TRUE){
      # Combine the specified columns from cData into one data.table
      benchGS <- rbind(cData[[3]][,c(1,4,5,7:13)],
                       cData[[4]][,c(1,4,5,7:13)],
                       cData[[3]][,c(1,4,5,14:20)],
                       cData[[4]][,c(1,4,5,14:20)])[, .SD[!duplicated(.SD)]]
      benchGS$mem <- as.numeric(benchGS$mem)
      
      # Add the new column rec_fold_rank and remove recordID, fold, and rank
      benchGS[, rec_fold_rank := paste0(recordID, "_", fold, "_", rank)]
      benchGS <- benchGS[, !c("recordID", "fold", "rank"), with = FALSE]
      bench <- rbind(benchOF, benchGS)
      
      # Susbet cData columns to remove benchmark data
      cData[["Genus"]] <- cData[["Genus"]][,1:6]
      cData[["Species"]] <- cData[["Species"]][,1:6]
      
      # Append bench data to cData
      cData <- list("cData"=cData, "benchData"=bench)
      rm(benchOF, benchGS, bench)
    } else {
      cData <- list("cData"=cData)
    }
  }
  
  # Post-Processing Strategy or MaxScore
  if(classStrategy == "Post-Processing" | classStrategy == "MaxScore"){
    # Iterate through all ranks in their totality (all models run)
    for (rank in subclassRank) {
      cData[[rank]] <- rbindlist(lapply(uniqFolds, function(foldD) {
        rbindlist(lapply(klenDat, function(klenD) {
          kTest <- if(klenD == "MixedKlen"){
            kTest <- kPosDataSub[(fold == foldD & klen != "1")]
          } else {
            kTest <- kPosDataSub[(fold == foldD & klen == klenD)]
          }
          modelNames <- names(models)[str_detect(names(models), rank)]
          .kphmm_classifyForwardCV(kTest, models[modelNames], uniqRec, rank, 
                                   foldD, subclassRank, klenD, dfTaxa, dbPath,
                                   classStrategy, runBenchmark, groupDict)
        }))
      }))
    }
    
    # If Post-processing, compute weighted average scores across ranks to determine best possible taxonomy
    if(classStrategy == "Post-Processing"){
      # Use weighted averaging of scores to determine best possible taxonomy across all levels
      cData <- .kphmm_weightAvgPost(recordDict, 
                                    cData[["Order"]], cData[["Family"]], 
                                    cData[["Genus"]], cData[["Species"]], 
                                    weights)
    }
    
    # If run benchmark is TRUE separate out benchmark data from classification data
    if(runBenchmark == TRUE){
      # Combine the specified columns from cData into one data.table
      bench <- rbind(cData[[1]][,c(1,4,5,7:13)],
                     cData[[2]][,c(1,4,5,7:13)],
                     cData[[3]][,c(1,4,5,7:13)],
                     cData[[4]][,c(1,4,5,7:13)],
                     cData[[1]][,c(1,4,5,14:20)],
                     cData[[2]][,c(1,4,5,14:20)],
                     cData[[3]][,c(1,4,5,14:20)],
                     cData[[4]][,c(1,4,5,14:20)])[, .SD[!duplicated(.SD)]]
      
      # Add the new column rec_fold_rank and remove recordID, fold, and rank
      bench[, rec_fold_rank := paste0(recordID, "_", fold, "_", rank)]
      bench <- bench[, !c("recordID", "fold", "rank"), with = FALSE]
      
      # Susbet cData columns to remove benchmark data
      cData[["Order"]] <- cData[["Order"]][,1:6]
      cData[["Family"]] <- cData[["Family"]][,1:6]
      cData[["Genus"]] <- cData[["Genus"]][,1:6]
      cData[["Species"]] <- cData[["Species"]][,1:6]
      
      # Append bench data to cData
      cData <- list("cData"=cData, "benchData"=bench)
      rm(bench)
    } else {
      cData <- list("cData"=cData)
    }
  }
  
  # Rbindlist cData across all ranks into one dataframe
  cData[["cData"]] <- rbindlist(cData[["cData"]])

  # If majority voting is TRUE
  if(majorityVoting == TRUE){
    # Seprate copy of original results
    cDataTmp <- data.frame(cData[["cData"]])
    
    # Create another data table for majority voting
    majVote <- cData[["cData"]]
    setDT(majVote)
    
    # Set keys for mvWt
    setkey(mvWt, "rank", "klen")
    
    # If not mixedklen, convert all klens >= 5 to ">=5"
    majVote <- majVote[klen != "MixedKlen",`:=`(klen = as.character(ifelse(as.numeric(klen) >= 5, ">=5", klen)))]
    setkey(majVote, "rank", "klen")
    
    # Merge weights with majVote
    majVote <- merge(mvWt, majVote)

    # Separate out z_scores and take a mean zscore grouping by rank, fold and group
    majVoteZ <- majVote[, .(z_score = mean(z_score)), by = .(fold, rank, group)]
    
    # Set keys for merge
    setkey(majVoteZ, "rank", "fold", "group")
    
    # Group by fold and rank and remove instance of no data, then sum up weighted votes
    majVote[, wtVote := ifelse(group == "No data", 0, weights * 1)] 
    majVote <- majVote[, .(voteSum = sum(wtVote)), by = .(fold, rank, group)][order(-voteSum)]
    majVote <- majVote[, .SD[1], by = .(fold, rank)] 
    
    # Set keys for merge
    setkey(majVote, "rank", "fold", "group")
    
    # Merge majVote and majVoteZ
    majVote <- merge(majVote, majVoteZ, by = c("rank", "fold", "group"))
    
    # Sort by z_score and choose the largest z_score group, then add MajorityVote and recordID
    majVote <- majVote[order(-z_score), .SD[1], by = .(fold, rank)]
    majVote[, `:=` (klen = "MajorityVote", recordID = uniqRec)]
    
    # rearrange columns and rbind back to cData
    cData[["cData"]] <- rbind(cDataTmp, majVote[,c("recordID","z_score","group","rank","fold","klen")]) %>% 
      as.data.table()
    
    # rm unneeded vars
    rm(majVote, majVoteZ, cDataTmp)
  }
  
  # Return cData
  return(cData)
}

# Internal function to find the previous taxonomic rank
.kphmm_prevRank <- function(rank) {
  switch(rank,
         "Family" = "Order",
         "Genus" = "Family",
         "Species" = "Genus")
}

# Internal function to look at the previously assigned taxa 
# and subset based on those previous classifications, if not, use models at rank above the higher rank
.kphmm_subsetbyPrevRank <- function(cData, models, recordDict, rank, uniqRec, foldD, klenD) {
  prevRankName <- .kphmm_prevRank(rank)
  prevGroups <- unique(cData[[prevRankName]][recordID == uniqRec & fold == foldD & klen == klenD]$group)
  if(prevGroups != "No data"){
    if(rank == "Family"){
      pattern <- paste0(rank, ";(", paste0(unique(recordDict[get(paste0(tolower(substr(prevRankName, 1, 3)), "ID")) %in% prevGroups]$famID), collapse = "|"), ")", sep="")
    }
    if(rank == "Genus"){
      pattern <- paste0(rank, ";(", paste0(unique(recordDict[get(paste0(tolower(substr(prevRankName, 1, 3)), "ID")) %in% prevGroups]$genID), collapse = "|"), ")", sep="")
    }
    if(rank == "Species"){
      pattern <- paste0(rank, ";(", paste0(unique(recordDict[get(paste0(tolower(substr(prevRankName, 1, 3)), "ID")) %in% prevGroups]$spID), collapse = "|"), ")", sep="")
    }
  } else {
    pattern <- "No data"
  }

  # If nothing is available, assign "No data", else subset model names by pattern above
  if(pattern=='' | pattern == "No data"){
    names <- "No data"
  } else {
    names <- names(models)[str_detect(names(models), pattern)]
  }
  return(names)
}

# Internal function to perform weighted averages to find best possible Order and Family - for Multi-Tiered Strategy only
.kphmm_weightAvgMulti <- function(recordDict, cDataOrd, cDataFam, weights){
  # Find all viable taxonomies that could exist from recordDict
  viableTaxonomies <- unique(paste(recordDict$ordID, "_", recordDict$famID, sep=""))
  viableTaxonomies <- str_split(viableTaxonomies, "_", simplify = TRUE)  # Convert to matrix
  
  # Rbindlist and filter by viable taxonomies from recordDict
  setkey(cDataOrd, "fold", "klen")
  setkey(cDataFam, "fold", "klen")
  
  # Performing a DT merge
  cDataOrdFam <- rbindlist(foreach(i=1:nrow(viableTaxonomies)) %do% {
    merge(cDataOrd[group == viableTaxonomies[i, 1]], cDataFam[group == viableTaxonomies[i, 2]], by = c("fold", "klen"), all = TRUE)
  })
  
  # Compute weighted average, retain all columns, and select the highest weighted average score
  cDataOrdFam[, weightedAvg := (z_score.x * weights[1] + z_score.y * weights[2]) / sum(weights)]
  cDataOrdFam <- cDataOrdFam[order(-weightedAvg), .SD[which.max(weightedAvg)], by = .(fold, klen)]
  
  # Return cData back to original format for genus and species runs
  cData <- list()
  cData[["Order"]] <- cDataOrdFam[,c(3, 4, 5, 6, 1, 2)]
  colnames(cData[["Order"]]) <- c("recordID","z_score","group","rank","fold","klen")
  cData[["Family"]] <- cDataOrdFam[,c(3, 8, 9, 10, 1, 2)]
  colnames(cData[["Family"]]) <- c("recordID","z_score","group","rank","fold","klen")
  
  # remove unneeded vars
  rm(cDataOrdFam, viableTaxonomies)
  
  # return cData
  return(cData)
}

# Internal function to perform weighted averages to find best possible taxonomy across all ranks - for Post Processing only
.kphmm_weightAvgPost <- function(recordDict, cDataOrd, cDataFam, cDataGen, cDataSp, weights){
  # Find all viable taxonomies that could exist from recordDict
  viableTaxonomies <- unique(paste(recordDict$ordID, ";", recordDict$famID, ";", recordDict$genID, ";", recordDict$spID, sep=""))
  viableTaxonomies <- str_split(viableTaxonomies, ";", simplify = TRUE)  # Convert to matrix
  
  # Rbindlist and filter by viable taxonomies from recordDict, then perform dt merges
  setkey(cDataOrd, "fold", "klen")
  setkey(cDataFam, "fold", "klen")
  setkey(cDataGen, "fold", "klen")
  setkey(cDataSp, "fold", "klen")
  
  # Perform data table merges across all ranks, filtering by viable taxonomies
  cDataAll <- suppressWarnings(rbindlist(foreach(i=1:nrow(viableTaxonomies)) %do% {
    Reduce(merge, list(cDataOrd[group == viableTaxonomies[i, 1]],
                       cDataFam[group == viableTaxonomies[i, 2]],
                       cDataGen[group == viableTaxonomies[i, 3]],
                       cDataSp[group == viableTaxonomies[i, 4]]))
  }))
  
  # Rename columns after merge
  colnames(cDataAll) <- c("fold", "klen", 
                          "recordID_1", "z_score_1", "group_1", "rank_1", 
                          "recordID_2", "z_score_2", "group_2", "rank_2", 
                          "recordID_3", "z_score_3", "group_3", "rank_3",
                          "recordID_4", "z_score_4", "group_4", "rank_4")
  
  # Compute weighted average, retain all columns, and select the highest weighted average score
  cDataAll[, weightedAvg := (z_score_1 * weights[1] + z_score_2 * weights[2] + z_score_3 * weights[3] + z_score_4 * weights[4]) / sum(weights)]
  cDataAll <- cDataAll[order(-weightedAvg), .SD[which.max(weightedAvg)], by = .(fold, klen)]
  
  # Return cData back to original format for genus and species runs
  cData <- list()
  cData[["Order"]] <- cDataAll[,c("recordID_1", "z_score_1", "group_1", "rank_1", "fold", "klen")]
  cData[["Family"]] <- cDataAll[,c("recordID_1", "z_score_2", "group_2", "rank_2", "fold", "klen")]
  cData[["Genus"]] <- cDataAll[,c("recordID_1", "z_score_3", "group_3", "rank_3", "fold", "klen")]
  cData[["Species"]] <- cDataAll[,c("recordID_1", "z_score_4", "group_4", "rank_4", "fold", "klen")]
  
  # Correct colnames
  foreach(i=1:length(cData)) %do% {
    colnames(cData[[i]]) <- c("recordID","z_score","group","rank","fold","klen")
  }
  
  # remove unneeded vars
  rm(viableTaxonomies, cDataAll)
  
  # return cData
  return(cData)
}

#kPosDataSub <- kPosData[(recordID == uniqRec[1])]
#models <- modelData$Models

#foldD <- "f1"
#klenD <- "MixedKlen"
#rank <- "Order"
#models <- models[str_detect(names(models), rank)]
#kPosDataSub <- if(klenD == "MixedKlen"){
#  kTest <- kPosDataSub[(fold == foldD & klen != "1")]
#} else {
#  kTest <- kPosDataSub[(fold == foldD & klen == klenD)]
#}

#.kphmm_classifyForwardCV(kPosDataSub, models, uniqRec, rank, foldD, subclassRank, klenD, dfTaxa, dbPath, classStrategy, runBenchmark, groupDict)

# Internal function to subset by position according to each model and use the forward algorithm on each model with each set of kmer data 
# (corresponding to one recordID each)
.kphmm_classifyForwardCV <- function(kPosDataSub, models, uniqRec, rank, foldD, subclassRank, klenD, dfTaxa, dbPath, classStrategy, runBenchmark, groupDict) {
  
  # Bench st
  benchS <- .bench_Func(runBenchmark, klenD, "kphmm_ClassifyCV", "None", 0, "Start")
  
  # Filter models by that specific klen and fold
  models <- Filter(Negate(is.null), lapply(models, function(m) m[[foldD]][[klenD]]))

  # If model length is greater than 1 - since one model doesnt tell you anything and has no point of comparison
  if(length(models)>1){
    # If not a MixedKlen model
    if(klenD != "MixedKlen"){
      # Actually run the test sequence on the models with the forward algorithm
      probsF <- unlist(foreach(j=1:length(models)) %do% {
        if(models[[j]]$modelData$success == TRUE){
          # Read in the model data
          model <- read_dbTable(dbPath, models[[j]]$modelData$model)
          modelPos <- model$alignment
          
          # group being compared
          groupDat <- models[[j]]$modelData$Group
          #print(paste0(foldD, "_", klenD, "_", groupDat, "_", j, sep=""))
            
          # If not running kmer length 1
          if(klenD != 1){
            # blastCheck db table
            blastCheck <- .kphmm_blastCheck(dfTaxa[(recordID == uniqRec), .(Id = recordID, Seq = dna)][1,], 
                                            groupDict[(group == groupDat & folds == foldD), .(Id = recordID, Seq = dna)])
              
            # if blastCheck does not have 0 rows, continue
            if(nrow(blastCheck)!=0){
              # Query match start and end
              qS <- as.numeric(blastCheck$QueryMatchStart)
              qE <- as.numeric(blastCheck$QueryMatchEnd)
                
              # trim by match start and end of train seqs (klen > 1 only)
              locations <- as.numeric(colnames(kPosDataSub[,2:(ncol(kPosDataSub)-2)]))
              locations <- locations[locations >= qS & locations <= qE]
                
              # Subset test sequence by those locations and renumber them from 1:ncol
              kTest <- kPosDataSub[, ..locations]
              colnames(kTest) <- as.character(1:ncol(kTest))
    
              # Finally subset by locations for that specific model
              locations <- colnames(modelPos)[colnames(modelPos) %in% colnames(kTest)]
              kTest <- kTest[, ..locations]
                
              # Run forward algorithm and normalize score by the length of the test seq (number of kmers in this case), name with group name for model
              suppressWarnings(setNames(forward(model, as.character(kTest[1,]), odds = TRUE)$score / length(as.character(kTest[1,])), groupDat))
                
            # otherwise forego the blast check, not ideal but at least still gives a score
            } else {
              # Else simply remove recordID, klen and fold columns
              locations <- as.numeric(colnames(kPosDataSub[,2:(ncol(kPosDataSub)-2)]))
              
              # Subset test sequence by those locations
              kTest <- kPosDataSub[, ..locations]
              
              # Subset by locations for that specific model
              locations <- colnames(modelPos)[colnames(modelPos) %in% colnames(kTest)]
              kTest <- kTest[, ..locations]
              
              # Run forward algorithm and normalize score by the length of the test seq (number of kmers in this case), name with group name for model
              suppressWarnings(setNames(forward(model, as.character(kTest[1,]), odds = TRUE)$score / length(as.character(kTest[1,])), groupDat))
            }
              
          # If running kmer length 1
          } else {
            # Else simply remove recordID, klen and fold columns
            locations <- as.numeric(colnames(kPosDataSub[,2:(ncol(kPosDataSub)-2)]))
              
            # Subset test sequence by those locations
            kTest <- kPosDataSub[, ..locations]
              
            # Run forward algorithm and normalize score by the length of the test seq (number of nucleotides in this case), name with group name for model
            suppressWarnings(setNames(forward(model, as.character(kTest[1,]), odds = TRUE)$score / length(as.character(kTest[1,])), groupDat))
          }
        # Else just set a score of -10
        } else { setNames(-10, models[[j]]$modelData$Group) }
      })
      
    # If a mixed length model
    } else {
      # Extract unique kmers from the model data for each unique group, then name with klen
      uniqModelK <- lapply(models, function(model) {
        kmers <- model$modelData$DiscrimKmers
        names(kmers) <- nchar(kmers)
        kmers
      })
      
      # Only proceed if there is data left after filtering by klen
      if(nrow(kPosDataSub) != 0) {
        kPosDataSub <- foreach(j = 1:length(uniqModelK)) %do% { 
          kPosDataSub[klen %in% unique(names(uniqModelK[[j]]))]
        }
      # Return an empty list if no data
      } else {
        kPosDataSub <- list()  
      }
      
      # Actually run the test sequence on the models with the forward algorithm
      probsF <- unlist(foreach(j=1:length(uniqModelK)) %do% {
        if(models[[j]]$modelData$success == TRUE){
            
          # group being compared
          groupDat <- models[[j]]$modelData$Group
          #print(paste0(foldD, "_", klenD, "_", groupDat, "_", j, sep=""))

          # blastCheck db table
          blastCheck <- .kphmm_blastCheck(dfTaxa[recordID == uniqRec, .(Id = recordID, Seq = dna)][1,], 
                                          groupDict[(group == groupDat & folds == foldD), .(Id = recordID, Seq = dna)])
            
          # if blastCheck does not have 0 rows, continue
          if(nrow(blastCheck)!=0){
            
            # Query match start and end
            qS <- as.numeric(blastCheck$QueryMatchStart)
            qE <- as.numeric(blastCheck$QueryMatchEnd)
              
            # trim by match start and end of train seqs (klen > 1 only)
            klenCol <- kPosDataSub[[j]][, "klen"]
            locations <- as.numeric(colnames(kPosDataSub[[j]][,2:(ncol(kPosDataSub[[j]])-2)]))
            locations <- locations[locations >= qS & locations <= qE]
              
            # Subset test sequence by those locations and renumber them from 1:ncol, add back a klength column
            kPosDataSub[[j]] <- kPosDataSub[[j]][, ..locations]
            colnames(kPosDataSub[[j]]) <- as.character(1:ncol(kPosDataSub[[j]]))
              
            # Split data table by klen and then cbind together
            kPosDataSub[[j]][, klen := klenCol]
            kPosDataSub[[j]] <- split(kPosDataSub[[j]], by = "klen")
            kPosDataSub[[j]] <- do.call(cbind, kPosDataSub[[j]])
              
            # Read in the model data
            model <- read_dbTable(dbPath, models[[j]]$modelData$model)
            modelPos <- model$alignment
              
            # Finally subset by kmer positions for that specific model
            colnames(modelPos) <- gsub("_",".",colnames(modelPos))
            locations <- colnames(modelPos)[colnames(modelPos) %in% colnames(kPosDataSub[[j]])]
            kPosDataSub[[j]] <- kPosDataSub[[j]][, ..locations]
              
            # Forward score normalized by num kmers in test sequence
            suppressWarnings(setNames(forward(model, as.character(kPosDataSub[[j]][1,]), odds = TRUE)$score / length(as.character(kPosDataSub[[j]][1,])), groupDat))
              
            # otherwise forego the blast check, not ideal but at least still gives a score
          } else {
            klenCol <- kPosDataSub[[j]][, "klen"]
            locations <- as.numeric(colnames(kPosDataSub[[j]][,2:(ncol(kPosDataSub[[j]])-2)]))
            kPosDataSub[[j]] <- kPosDataSub[[j]][, ..locations]
            
            # Split data table by klen and then cbind together
            kPosDataSub[[j]][, klen := klenCol]
            kPosDataSub[[j]] <- split(kPosDataSub[[j]], by = "klen")
            kPosDataSub[[j]] <- do.call(cbind, kPosDataSub[[j]])
            
            # Read in the model data
            model <- read_dbTable(dbPath, models[[j]]$modelData$model)
            modelPos <- model$alignment
            
            # Finally subset by kmer positions for that specific model
            colnames(modelPos) <- gsub("_",".",colnames(modelPos))
            locations <- colnames(modelPos)[colnames(modelPos) %in% colnames(kPosDataSub[[j]])]
            kPosDataSub[[j]] <- kPosDataSub[[j]][, ..locations]
            
            # Forward score normalized by num kmers in test sequence
            suppressWarnings(setNames(forward(model, as.character(kPosDataSub[[j]][1,]), odds = TRUE)$score / length(as.character(kPosDataSub[[j]][1,])), groupDat))
          }
        # Else just set a score of -10
        } else { setNames(-10, models[[j]]$modelData$Group) }
      })
    }
    
    # Keep forward scores > -10
    probsF <- probsF[probsF>-10]
    
    # If length probsF is not 0
    if(length(probsF)!=0){
      # Separate out positive normalized infinite forward scores, remove negative normalized infinite scores
      # If +inf present, assign a z_score 1 higher than highest numerical z_score to that group since this would represent the highest possible score 
      z_scoreInf <- probsF[is.infinite(probsF) & probsF > 0]
      z_scoreInf[is.infinite(z_scoreInf)] <- max(probsF) + 0.1
      probsF <- probsF[!is.infinite(probsF)]
      
      # remove unneeded vars
      suppressWarnings(rm(kTestSubLoc, modelPos, model, uniqueLoc))
      
      # Determine z-scores from the normalized forward scores
      # Calculate mean and standard deviation
      mean_score <- mean(probsF)
      sd_score <- sd(probsF)
      
      # Calculate z-scores
      z_scores <- (probsF - mean_score) / sd_score
      
      # If z_scoreinf is not empty, append to z-scores calculated
      if(length(z_scoreInf) != 0){
        z_scores <- append(z_scoreInf, z_scores)
      }
      
      # rm unneeded vars
      suppressWarnings(rm(kTestSub, mean_score, sd_score))
      
      # Convert results data to dataframe format, rbindlist and then
      # select the top scoring group based on odds score from the forward algorithm
      kTestProb <- rbindlist(foreach(j=1:length(z_scores)) %do% { 
        data.table("recordID" = uniqRec,
                   "z_score" = as.numeric(z_scores[j]),
                   "group" = names(z_scores)[j], 
                   "rank" = rank, 
                   "fold" = foldD, 
                   "klen" = klenD)
      })
      
      # Sort by highest z-score and pick the top one (depending on classification strategy)
      if(classStrategy == "MaxScore" | classStrategy == "Hierarchical"){
        kTestProb <- kTestProb[order(-z_score)][1]
      }
      if(classStrategy == "Multi-Tiered"){
        if("Genus" %in% subclassRank | "Species" %in% subclassRank){
          kTestProb <- kTestProb[order(-z_score)][1]
        }
      }
      
    # Else simply return a dataframe with no data for group and -10 for z-score
    } else {
      kTestProb <- data.table("recordID" = uniqRec,
                              "z_score" = -10,
                              "group" = "No data", 
                              "rank" = rank, 
                              "fold" = foldD, 
                              "klen" = klenD)
    }
    
  # Else simply return a dataframe with no data for group and -10 for z-score
  } else {
    kTestProb <- data.table("recordID" = uniqRec,
                            "z_score" = -10,
                            "group" = "No data", 
                            "rank" = rank, 
                            "fold" = foldD, 
                            "klen" = klenD)
  }

  # Bench end
  benchE <- .bench_Func(runBenchmark, klenD, "kphmm_ClassifyCV", "None", 0, "End")
  
  # If there is benchmarking data
  if(runBenchmark == TRUE){
    kTestProb <- cbind(kTestProb, benchS, benchE)
  }
  
  # Return kTestProb dataframe
  return(kTestProb)
}

# Internal function to check with BLAST if start and end location sites differ between test seq and train seqs
# and if they do trim test sequence by the query match start and end of the top ranked hit
.kphmm_blastCheck <- function(recSeq, groupDict){
  # BLAST all test records against all sequences belonging to a particular group, set minIdentity a bit lower to accomodate dissimilar sequences
  blastResults <- blaster::blast(recSeq, groupDict, maxAccepts = 1, minIdentity = 0.6, output_to_file = FALSE)
  # Select relevant columns
  blastResults <- blastResults[,c(1,3:4)]
  colnames(blastResults)[1] <- "recordID"
  blastResults$QueryMatchStart <- as.numeric(blastResults$QueryMatchStart)
  blastResults$QueryMatchEnd <- as.numeric(blastResults$QueryMatchEnd)
  return(blastResults)
}

# Main function to classify sequences independently from sequences provided by the user on trained kphmm models (which may or may have not been trained by the user)
# modelResults <- kphmm_ClassifyI(subclassRank){
# }

##### Benchmarking Functions #####

# Internal function to benchmark time and memory for functions
.bench_Func <- function(runBenchmark, klen, funcName, subFunc, dbSize, stORend){
  if(runBenchmark == TRUE){
    # Pryr for memory used prior to function run, system.time for time start, file.info for database size
    mem <- as.numeric(pryr::mem_used())
    time <- Sys.time()
    # Combine in a DT and return
    bench <- data.table("mem"=mem, "time"=time, "dbSize"=dbSize, "klen"=klen, "funcName"=funcName, "subFuncName"=subFunc, "stORend"=stORend)
  } else {
    bench <- NA
  }
  return(bench)
}

# Internal function to aggregate benchmark results
.bench_Aggregate <- function(runBenchmark, bench, nSeqs, numCores, numFolds, func){
  if(runBenchmark == TRUE){
    # If running function kphmm_ClassifyCv
    if(func == "kphmm_ClassifyCV"){
      # Perform the necessary operations on the data.table
      bench <- bench[, .(`memUsage (Mb)` = round(as.numeric((mem[stORend == 'End'] - mem[stORend == 'Start']) / 1024^2), 8),
                         `elapsedTime (mins)` = round(as.numeric(difftime(time[stORend == 'End'], time[stORend == 'Start'], units = "mins")), 8),
                         `dbUsage (Mb)` = round(as.numeric(dbSize[stORend == 'End'] - dbSize[stORend == 'Start']), 8)),
                     by = .(klen, rec_fold_rank)
      ]
      
      # Replace negative memory usage with 0
      bench[`memUsage (Mb)` < 0, `memUsage (Mb)` := 0]
      
      # Aggregate by 'klen' - take total time and memory
      bench <- bench[, .(`memUsage Total (Mb)` = sum(`memUsage (Mb)`),
                         `elapsedTime Total (mins)` = sum(`elapsedTime (mins)`)),
                     by = .(klen)
      ]
      
      # Take time and memory per seq
      bench <- bench[, .(`memUsage AvgPerSeq (Mb)` = (`memUsage Total (Mb)` / nSeqs),
                         `elapsedTime AvgPerSeq (mins)` = `elapsedTime Total (mins)` / nSeqs),
                     by = .(klen, `memUsage Total (Mb)`,`elapsedTime Total (mins)`)
      ]
      # Add back funcName
      bench[, funcName := func]
      
    # If running function kphmm_Train
    } else if(func == "kphmm_Train"){
      # Convert the list of data.tables to a data.table
      bench <- rbindlist(bench, idcol="groupNum")
      
      # Perform the necessary operations on the data.table
      bench <- bench[, .(`memUsage (Mb)` = round(as.numeric((mem[stORend == 'End'] - mem[stORend == 'Start']) / 1024^2), 8),
                         `elapsedTime (mins)` = round(as.numeric(difftime(time[stORend == 'End'], time[stORend == 'Start'], units = "mins")), 8),
                         `dbUsage (Mb)` = round(as.numeric(dbSize[stORend == 'End'] - dbSize[stORend == 'Start']), 8)),
                     by = .(klen, funcName, subFuncName, groupNum)
      ]
      
      # Replace negative memory usage with 0
      bench[`memUsage (Mb)` < 0, `memUsage (Mb)` := 0]
      
      # count number of observations for each klen
      numObs <- c(table(bench$klen))
      
      # Aggregate by 'klen' and sum metrics and divide by number of unique klengths
      bench <- bench[, .(`memUsage (Mb)` = sum(`memUsage (Mb)`) / length(numObs),
                         `elapsedTime (mins)` = sum(`elapsedTime (mins)`) / length(numObs),
                         `dbUsage (Mb)` = sum(`dbUsage (Mb)`)),
                     by = .(klen, funcName, subFuncName)]
      # Split by klen
      bench <- split(bench, by=c("klen"))
      
      # Sort numObs by bench names
      numObs <- numObs[order(match(names(numObs),names(bench)))]
      
      # Divide values by number of observations for each klen and multiply by mean number of observations 
      # to get standardized values
      foreach(i=1:length(numObs)) %do% {
        bench[[i]] <- bench[[i]][, .(`memUsage (Mb)` = (`memUsage (Mb)` / numObs[i]) * mean(numObs),
                                     `elapsedTime (mins)` = (`elapsedTime (mins)` / numObs[i]) * mean(numObs),
                                     `dbUsage (Mb)` = (`dbUsage (Mb)` / numObs[i]) * mean(numObs)),
                                 by = .(klen, funcName, subFuncName)]
      }
      
      # Rbindlist again and add back func names columns
      bench <- rbindlist(bench)
    
    # If using other function besides kPos_Mat
    } else {
      # Convert the list of data.tables to a data.table
      bench <- rbindlist(bench)
      
      # Perform the necessary operations on the data.table
      bench <- bench[, .(`memUsage (Mb)` = round(as.numeric((mem[stORend == 'End'] - mem[stORend == 'Start']) / 1024^2), 8),
                         `elapsedTime (mins)` = round(as.numeric(difftime(time[stORend == 'End'], time[stORend == 'Start'], units = "mins")), 8),
                         `dbUsage (Mb)` = round(as.numeric(dbSize[stORend == 'End'] - dbSize[stORend == 'Start']), 8)),
                     by = .(klen, funcName, subFuncName)
      ]
      
      # Replace negative memory usage with 0
      bench[`memUsage (Mb)` < 0, `memUsage (Mb)` := 0]
      
      # Aggregate by 'klen'
      bench <- bench[, .(`memUsage (Mb)` = sum(`memUsage (Mb)`),
                         `elapsedTime (mins)` = sum(`elapsedTime (mins)`),
                         `dbUsage (Mb)` = sum(`dbUsage (Mb)`)),
                     by = .(klen, funcName, subFuncName)
      ]
    }
    
    # Add additional columns with system information
    bench[, nSeqsUsed := nSeqs]
    bench[, nFoldsUsed := numFolds]
    bench[, numberCoresUsed := numCores]
    bench[, platform := benchmarkme::get_platform_info()$OS.type[1]]
    bench[, arch := benchmarkme::get_platform_info()$r_arch[1]]
    bench[, rVersion := benchmarkme::get_r_version()$version.string[1]]
    bench[, totalRam := capture.output(print(benchmarkme::get_ram()))]
    bench[, cpu := benchmarkme::get_cpu()$model_name[1]]
    
    # Remove duplicate rows if any
    bench <- unique(bench)
    
  # Else simply assign NA  
  } else {
    bench <- NA
  }
  return(bench)
}

####### Results Plotting Function ########

# Main function for writing out results plots to png files
plot_Results <- function(kmerData, kmerResults, modelData, modelResults, 
                         klen, numFolds, classRank, subclassRank, 
                         runBenchmark, numCores, dpi){
  
  # Create a directory for results plots
  dirP <- paste0(dir, "_Plots")
  dir.create(dirP, showWarnings = FALSE)
  
  # All klens appended with MixedKlen
  klenDat <- append(klen, "MixedKlen")
  
  # List for filepaths of all plots generated
  plotNames <- list()
  
  # Message update
  cat(green(paste0("Writing out accuracy curves to ", dirP, "...\n", sep="")))
  
  ### ROC curves ###
  
  # rbindlist across klens and rank for AUC values
  auc <- rbindlist(foreach(i=1:length(modelResults$AccMetrics$`ROC-AUC`)) %do% {
    rbindlist(foreach(j=1:length(modelResults$AccMetrics$`ROC-AUC`[[i]])) %do% {
      data.frame("auc" = modelResults$AccMetrics$`ROC-AUC`[[i]][[j]],
                 "klen" = names(modelResults$AccMetrics$`ROC-AUC`[[i]])[j],
                 "rank" = names(modelResults$AccMetrics$`ROC-AUC`)[i])
    })
  })
  
  # rbindlist across klens and rank for the ROC data
  roc <- rbindlist(foreach(i=1:length(subclassRank)) %do% {
    rbindlist(modelResults$AccMetrics$Data[[i]], idcol="klen") %>% 
      mutate(rank = subclassRank[i])
  })
  
  # Scale x.sorted values from 0-1 to 1 when grouped by klen and rank
  roc <- roc %>%
    # drop NaN or Inf values
    filter(x.sorted != "Inf") %>%
    # group by fold and klen
    group_by(klen, rank) %>%
    # For the other curves (besides ROC) scale confidence values 0-1
    mutate(scaledXsorted = .scaleValues(x.sorted))
  
  # setDTs, keys and merge
  setDT(roc)
  setDT(auc)
  setkey(roc, "klen", "rank")
  setkey(auc, "klen", "rank")
  roc <- merge.data.table(roc, auc) %>% mutate(klen_auc = paste0(klen, ", ", auc, sep=""))
  
  # Plot out ROC curves per rank color-coding by klen
  ggROC <- list()
  foreach(i=1:length(subclassRank)) %do% {
    ggROC[[i]] <- suppressWarnings(ggplot(data=roc %>% filter(rank == subclassRank[i]), aes(x=fpr, y = tpr, col = klen_auc)) + 
                                     geom_line(size=1.2, alpha=0.8) + 
                                     scale_colour_viridis_d(name = "kLength & AUC") + 
                                     xlab("Average FPR") + ylab("Average TPR")  + 
                                     ggtitle(paste0("Class: ", classRank, 
                                                    ", Rank: ", subclassRank[i],
                                                    ", nFolds = ", numFolds, 
                                                    ", nGroups = ", length(unique(modelResults$`CorrectID%`[[1]]$actualTaxa)), 
                                                    ", nSeqs = ", length(unique(modelResults$CorrectIDs[[1]]$recordID)), sep="")) +
                                     theme(text = element_text(size=16)) +
                                     annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + 
                                     theme(panel.spacing = unit(1, "lines")))
  }
  
  ### MCC curves ###
  
  # Remove any NAs for MCC plots, and calaculate means over intervals of 0.02 scaled score
  MCC <- roc %>%
    filter(!is.na(MCC_score)) %>% 
    group_by(klen, rank) %>%
    mutate(int = cut_interval(scaledXsorted, 
                              n = (max(scaledXsorted)-min(scaledXsorted))/0.02)) %>%
    group_by(klen, rank, int) %>%
    dplyr::summarise(mean(MCC_score)) %>%
    separate(int, c("St", "End"), sep=",") %>% 
    ungroup() %>% 
    group_by(klen, rank)
  MCC$int <- as.numeric(gsub("]", "", MCC$End))
  
  # Plot out MCC curves per rank color-coding by klen
  ggMCC <- list()
  foreach(i=1:length(subclassRank)) %do% {
    ggMCC[[i]] <- suppressWarnings(ggplot(data=MCC %>% filter(rank == subclassRank[i]), aes(x=int, y = `mean(MCC_score)`, col = klen)) + 
                                     geom_line(size=1.2, alpha=0.8) + 
                                     scale_colour_viridis_d(name = "kLength") + 
                                     xlab("Confidence Threshold") + ylab("Average MCC score")  +
                                     theme(text = element_text(size=16)) +
                                     theme(panel.spacing = unit(1, "lines")))
  }
  
  ### F1 curves ###
  
  # Remove any NAs for F1 plots
  F1 <- roc %>%
    filter(!is.na(F1_score))
  
  ggF1 <- list()
  foreach(i=1:length(subclassRank)) %do% {
    ggF1[[i]] <- suppressWarnings(ggplot(data=F1 %>% filter(rank == subclassRank[i]), aes(x=scaledXsorted, y = F1_score, col = klen)) + 
                                    geom_line(size=1.2, alpha=0.8) + 
                                    scale_colour_viridis_d(name = "kLength") + 
                                    xlab("Confidence Threshold") + ylab("Average F1 score")  +
                                    theme(text = element_text(size=16)) +
                                    theme(panel.spacing = unit(1, "lines")))
  }
  
  ### Cohens-Kappa curves ###
  
  # Remove any NAs for CK plots
  CK <- roc %>%
    filter(!is.na(F1_score))
  
  ggCK <- list()
  foreach(i=1:length(subclassRank)) %do% {
    ggCK[[i]] <- suppressWarnings(ggplot(data=CK %>% filter(rank == subclassRank[i]), aes(x=scaledXsorted, y = cohens_kappa, col = klen)) + 
                                    geom_line(size=1.2, alpha=0.8) + 
                                    scale_colour_viridis_d(name = "kLength") + 
                                    xlab("Confidence Threshold") + ylab("Average cohens-kappa score")  +
                                    theme(text = element_text(size=16)) +
                                    theme(panel.spacing = unit(1, "lines")))
  }
  
  ### Combined ROC/MCC/CK plot ###
  
  # Grid arrangement of all plots (ROC/MCC/F1/CK) combined
  ggComb <- list()
  foreach(i=1:length(subclassRank)) %do% {
    ggComb[[i]] <- grid.arrange(ggROC[[i]], ggMCC[[i]], ggF1[[i]], ggCK[[i]], nrow = 4, ncol = 1)
    ggsave(paste0(dirP, "/Acc_Curves_", subclassRank[i], ".png", sep=""), 
           plot=ggComb[[i]], 
           width = 10, height = 15, dpi = dpi, bg="white")
  }
  
  # File names for accuracy curve plots
  plotNames[["AccCurves"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "Acc_Curves")], sep="")
  
  # rm unneeded vars
  rm(auc, roc, ggComb, CK, ggCK, F1, ggF1, MCC, ggMCC, ggROC)
  
  ### Color-coded phylogenetic trees of percentage correct ids ###
  
  # if running at all ranks
  if(length(subclassRank) == 4){
    
    # Message update
    cat(green(paste0("Writing out %CorrectId plots to ", dirP, "...\n", sep="")))
    
    # rbind results from modelResults
    taxa <- rbindlist(modelResults$`CorrectID%`, idcol="rank")
    taxa <- taxa
    colnames(taxa)[2] <- "group"
    taxa <- taxa %>% group_by(rank) %>% group_split()
    taxaChar <- c("famID", "genID", "ordID", "spID")
    
    # Perform the operations using data.table syntax within the foreach loop
    taxa <- foreach(i=1:length(taxa)) %do% {
      # Determine the merge column and parent column based on the value of i
      merge_col <- switch(i,
                          "1" = "famID",
                          "2" = "genID",
                          "3" = "ordID",
                          "spID")
      
      # Merge with dfTaxa and rename/reorganize columns
      merge(taxa[[i]], dfTaxa, by.x = "group", by.y = merge_col) %>%
        mutate(actGroup = group) 
    }
    
    # Assign colnames with taxaChar and select the right columns
    foreach(i=1:length(taxa)) %do% {
      colnames(taxa[[i]])[1] <- taxaChar[i]
      taxa[[i]] <- taxa[[i]] %>% select(c(actGroup, `CorrectId%`, rank, fold, klen, n, spID, genID, famID, ordID, classID))
    }
    
    # Rbindlist taxa
    taxa <- rbindlist(taxa)
    
    # For class percentage correct ids per klength
    classPercent <- taxa %>% 
      select(c(`CorrectId%`, klen, ordID, n)) %>%
      group_by(klen) %>%
      mutate(trait = mean(`CorrectId%`)) %>%
      select(!c(`CorrectId%`, ordID)) %>%
      dplyr::distinct() %>%
      mutate(nSum = sum(n)) %>%
      select(!c(n)) %>%
      dplyr::distinct() %>%
      group_split()
    
    # Make a list per klength
    classPercentL <- list()
    foreach(i=1:length(klenDat)) %do% {
      classPercentL[[i]] <- tibble(label = "Mammalia",
                                   trait = classPercent[[i]]$trait[1],
                                   nSum = classPercent[[i]]$nSum[1])
    }
    
    # List for phylogenetic tree plots
    ggTreePlot <- replicate(length(klenDat), list())
    
    # For each klength make a phylogenetic tree of % correct ids
    foreach(i=1:length(klenDat)) %do% {
      
      # Filter by fold and klen and order and calculate correct ID percentages
      taxatest <- taxa %>% 
        filter(klen == klenDat[i]) %>%
        select(c(actGroup, fold, n, `CorrectId%`, classID, ordID, famID, genID, spID)) %>% 
        mutate(label = actGroup) %>%
        select(!actGroup) %>%
        dplyr::distinct() 
      
      # For summing number of classifications per taxonomic group
      taxatestn <- taxatest %>%
        select(c(label, n, fold)) %>%
        dplyr::distinct() %>%
        group_by(label) %>%
        dplyr::summarise(nSum = sum(n)) %>%
        ungroup() %>%
        as.data.table()
      setkey(taxatestn, label)
      
      # For taking the mean of correct id percentages across folds per group
      taxatestmean <- taxatest %>%
        select(c(label, `CorrectId%`, fold)) %>%
        dplyr::distinct() %>%
        group_by(label) %>%
        dplyr::summarise(trait = mean(`CorrectId%`)) %>%
        ungroup() %>%
        as.data.table()
      setkey(taxatestmean, label)
      
      # Merge them all together
      taxatest <- taxatest %>% 
        select(c(label, classID, ordID, famID, genID, spID)) %>%
        as.data.table()
      setkey(taxatest, label)
      taxatest <- Reduce(merge, list(taxatest,
                                     taxatestn,
                                     taxatestmean)) %>% dplyr::distinct()
      
      # rm unneeded vars
      rm(taxatestmean, taxatestn)
      
      # Split up larger orders based on size (1000 obs or higher)
      tallyGroups <- taxatest %>% 
        group_by(ordID) %>% 
        tally()
      ordGroups <- tallyGroups %>% 
        group_by(n >= 1000) %>% 
        group_split()
      
      # If groups exist larger than 100 obs
      if(length(ordGroups)>1){
        ordSmall <- ordGroups[[1]] %>% group_by(1:n()) %>% group_split()
        ordSmall <- list(unlist(foreach(l=1:length(ordSmall)) %do% ordSmall[[l]]$ordID))
        ordLarge <- ordGroups[[2]] %>% group_by(1:n()) %>% group_split()
        ordLarge <- foreach(l=1:length(ordLarge)) %do% ordLarge[[l]]$ordID
        ordGroupsList <- append(ordLarge, ordSmall)
        # If only smaller groups exist
      } else {
        ordGroupsList <- list(as.character(tallyGroups$ordID))
      }
      
      # Sublists for each ggTree plot
      ggTreePlot[[i]] <- replicate(length(ordGroupsList), list())
      
      # Rm unneeded vars
      suppressWarnings(rm(ordSmall, ordLarge, ordGroups, tallyGroups))
      
      # For each order group
      foreach(j=1:length(ordGroupsList)) %do% {
        
        # Filter by order being used and remove taxonomy columns
        taxatest2 <- taxatest %>%
          filter(ordID %in% ordGroupsList[[j]]) %>%
          select(!c(classID, ordID, famID, genID, spID)) %>%
          dplyr::distinct()
        taxatest2 <- as_tibble(rbind(taxatest2, classPercentL[[i]]))
        
        # Create additional dataframes to determine which families, orders, 
        # species and genuses will be included in the tree
        taxaSubSp <- taxa %>% 
          select(!c(actGroup, fold, `CorrectId%`, n, rank)) %>%
          filter(spID %in% taxatest2$label & klen == klenDat[i] & ordID %in% ordGroupsList[[j]])
        taxaSubGen <- taxa %>% 
          select(!c(actGroup, fold, `CorrectId%`, n, rank)) %>%
          filter(genID %in% taxatest2$label & klen == klenDat[i] & ordID %in% ordGroupsList[[j]])
        
        # If some genuses were not found in sp dataframe add those to the sp dataframe
        if(nrow(taxaSubGen %>% filter(!taxaSubGen$genID %in% taxaSubSp$genID))>0){
          taxaSub <- rbind(taxaSubGen %>% filter(!taxaSubGen$genID %in% taxaSubSp$genID), taxaSubSp)
          # Place "N/A" for species when genuses posses no species level trait data
          taxaSub$spID <- ifelse(!taxaSub$spID %in% taxaSubSp$spID, "N/A", taxaSub$spID)
        } else {
          taxaSub <- taxaSubSp
        }
        
        # Rm unneeded vars
        rm(taxaSubSp, taxaSubGen)
        
        # Remove unneeded cols from taxaSub and reorganize cols
        taxaSub <- taxaSub %>% 
          select(!c(klen)) %>% 
          dplyr::distinct()
        taxaSub <- taxaSub[,c(5,4,3,2,1)]
        
        # Convert to tree object with taxonomizr, first filtering by species and genuses
        colnames(taxaSub)[1:5] <- c("class", "order", "family", "genus", "species")
        taxaSub <- as.matrix(taxaSub)
        tree <- makeNewick(taxaSub)
        
        # Convert to phylo object
        taxaPhylo <- ape::read.tree(text = tree)
        
        # Convert to tibble 
        taxaTibb <- as_tibble(taxaPhylo)
        
        # Merge tree data with trait data
        taxaTibbMerge <- full_join(taxaTibb, taxatest2, by = c('label')) %>%
          mutate(label = ifelse(label!="N/A" & label != "", 
                                ifelse(!is.na(nSum),
                                       gsub("_", " ", paste0(label, " (n=", nSum, ")", sep="")),
                                       gsub("_", " ", paste0(label, " (n=", 0, ")", sep=""))),
                                label)) %>%
          mutate(nSum = ifelse(is.na(nSum),
                               1,
                               as.numeric(nSum))) %>% as.treedata()
        
        # ggTree plot generate per klength
        ggTreePlot[[i]][[j]] <- ggtree(taxaTibbMerge, layout = "circular") + 
          geom_nodelab(geom = 'label', aes(fill = trait, size=nSum), fontface = "bold") +
          geom_tiplab(geom = 'label', aes(fill = trait, size=nSum), fontface = "bold") +
          scale_fill_continuous_sequential(palette="Red-Blue", rev=TRUE, p1=1, p2=1, alpha=0.8, name = "Avg Correct Id%") +
          theme_tree("gray80") + 
          scale_size_continuous(range = c(6.5, 11.5)) +  
          guides(size = "none")
        
        # Tree rotation - suppressMessages for annoying console outputs
        suppressMessages(ggTreePlot[[i]][[j]] <- rotate_tree(ggTreePlot[[i]][[j]], -40))
        
        # If more than one plot grouping
        if(length(ordGroupsList)>1){
          ggsave(paste0(dirP, "/CorrectId_Tree", "_klen", klenDat[i], "_plot", j, ".png", sep=""), 
                 plot=ggTreePlot[[i]][[j]] +
                   theme(legend.position = 'top',
                         plot.title = element_text(size = 55, hjust=0.5),
                         legend.key.size = unit(4, 'cm'),
                         legend.title=element_text(size=40),
                         legend.text=element_text(size=35)) +
                   ggtitle(paste0("Order(s): ", paste0(ordGroupsList[[j]], collapse=", "), ", klength: ", klenDat[i],  ", nFolds = ", numFolds, sep="")), 
                 width = 49, height = 49, dpi = dpi, bg="white")
          
          # If only one plot grouping
        } else {
          ggsave(paste0(dirP, "/CorrectId_Tree", "_klen", klenDat[i], "_plot", j, ".png", sep=""), 
                 plot=ggTreePlot[[i]][[j]] +
                   theme(legend.position = 'top',
                         plot.title = element_text(size = 55, hjust=0.5),
                         legend.key.size = unit(4, 'cm'),
                         legend.title=element_text(size=40),
                         legend.text=element_text(size=35)) +
                   ggtitle(paste0("Class: ", classRank, ", klength: ", klenDat[i], ", nFolds = ", numFolds, sep="")), 
                 width = 49, height = 49, dpi = dpi, bg="white")
        }
        # Rm unneeded vars
        rm(tree, taxaPhylo, taxaTibb, taxaTibbMerge, taxatest2)
      }
    }
    
    # File names for accuracy curve plots
    plotNames[["CorrectId%"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "CorrectId_Tree")], sep="")
    
    # Trim whitespace in the margins of trees 
    foreach(i=1:length(plotNames[["CorrectId%"]])) %do% {
      image_write(image_trim(image_read(plotNames[[2]][i])), path = plotNames[[2]][i], format = "png")
      gc()
    }
    
    # Rm unneeded vars
    rm(ggTreePlot, ordGroupsList, classPercent, classPercentL, taxatest, taxa)
  }
  
  ### Bar plots of ranks/folds/klens lacking discriminative kmers (as a percentage of total number of groups) ###
  
  # Message update
  cat(green(paste0("Writing out barplots to ", dirP, "...\n", sep="")))
  
  # If groups were found to possess no discrim kmers
  if(nrow(kmerResults[["kFail"]])>0){
    # Unique folds present in kFail
    uniqFolds <- unique(kmerResults[["kFail"]]$fold)
    
    # Tally the number of groups per rank/fold/klen with no discrim kmers
    kFail <- kmerResults[["kFail"]] %>% 
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nFail = n()) %>%
      as.data.table()
    setkey(kFail, "rank", "fold", "klen")
    
    # Then grab the kTally results for groups that do have discrim kmers
    kSuccess <- kmerResults[["kTally"]] %>% 
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nSuccess = n()) %>%
      as.data.table()
    setkey(kSuccess, "rank", "fold", "klen")
    
    # Merge kFail with kSuccess
    kFS <- merge(kFail, kSuccess) %>%
      rowwise() %>%
      mutate(`kFail%` = (nFail / (nFail + nSuccess)) * 100) %>%
      mutate(nSum = nFail + nSuccess) %>%
      mutate(`Rank & Fold` = paste0(rank, " ", fold, sep=""))
    
    # Setting up y-axis label names for plot
    yLabels <- setNames(kFS$`Rank & Fold`, kFS$rank)
    
    # Set the labels as a factor for proper ordering by rank
    kFS <- kFS %>% 
      mutate(`Rank & Fold` = factor(`Rank & Fold`, 
                                    levels = c(unique(yLabels[names(yLabels)=="Species"]),
                                               unique(yLabels[names(yLabels)=="Genus"]), 
                                               unique(yLabels[names(yLabels)=="Family"]), 
                                               unique(yLabels[names(yLabels)=="Order"])))) %>%
      as.data.table()
    
    # Set up plot title
    pTitle <- paste0("Class: ", classRank, sep="")
    if(any(str_detect(kFS$rank, "Order"))){
      pTitle <- paste0(pTitle, ", nOrd=", ceiling(mean(kFS[(rank=="Order")]$nSum)), sep="")
    }
    if(any(str_detect(kFS$rank, "Family"))){
      pTitle <- paste0(pTitle, ", nFam=", ceiling(mean(kFS[(rank=="Family")]$nSum)), sep="")
    }
    if(any(str_detect(kFS$rank, "Genus"))){
      pTitle <- paste0(pTitle, ", nGen=", ceiling(mean(kFS[(rank=="Genus")]$nSum)), sep="")
    }
    if(any(str_detect(kFS$rank, "Species"))){
      pTitle <- paste0(pTitle, ", nSp=", ceiling(mean(kFS[(rank=="Species")]$nSum)), sep="")
    }
    
    # Create a ggplot bar plot
    ggkFail <- ggplot(kFS, aes(x = as.ordered(`Rank & Fold`), y = `kFail%`, 
                               fill = as.factor(klen))) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "Rank & Fold", 
           y = "Percentage of taxonomic groups without discrim kmers") +
      coord_flip() +
      theme_minimal() +
      theme(plot.title = element_text(size = 25, hjust=0.5),
            legend.position = 'top',
            legend.key.size = unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=18),
            axis.title = element_text(size = (22)),
            axis.text = element_text(size = (18))) +
      ggtitle(pTitle)
    
    # Write out to a file in the dirP directory
    ggsave(paste0(dirP, "/Barplot_kFail.png", sep=""), 
           plot=ggkFail, 
           width = 12, height = 17, dpi = dpi, bg="white")
    
    # File names for bar plot
    plotNames[["kmerFail%"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "kFail")], sep="")
    
    # Rm unneeded vars
    rm(kFail, kSuccess, kFS, yLabels, ggkFail)
  }
  
  ### Bar plot of groups (per fold/rank/klen) failing during training (with discrim kmers) ###
  
  # If groups were found to fail during training
  if(nrow(modelData[["Fail"]])>0){
    # Unique folds present in training fail
    uniqFolds <- unique(modelData[["Fail"]]$fold)
    
    # Tally the number of groups per rank/fold/klen where training failed
    tFail <- modelData[["Fail"]] %>% 
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nFail = n()) %>%
      as.data.table()
    setkey(tFail, "rank", "fold", "klen")
    
    # Then grab the kTally results for groups that do have discrim kmers and thus had models trained
    kSuccess <- kmerResults[["kTally"]] %>% 
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nSuccess = n()) %>%
      mutate(klen = as.character(klen)) %>%
      as.data.table()
    setkey(kSuccess, "rank", "fold", "klen")
    
    # Merge kFail with kSuccess
    tFS <- merge(tFail, kSuccess) %>%
      rowwise() %>%
      mutate(`tFail%` = (nFail / (nFail + nSuccess)) * 100) %>%
      mutate(nSum = nFail + nSuccess) %>%
      mutate(`Rank & Fold` = paste0(rank, " ", fold, sep=""))
    
    # Setting up y-axis label names for plot
    yLabels <- setNames(tFS$`Rank & Fold`, tFS$rank)
    
    # Set the labels as a factor for proper ordering by rank
    tFS <- tFS %>% 
      mutate(`Rank & Fold` = factor(`Rank & Fold`, 
                                    levels = c(unique(yLabels[names(yLabels)=="Species"]),
                                               unique(yLabels[names(yLabels)=="Genus"]), 
                                               unique(yLabels[names(yLabels)=="Family"]), 
                                               unique(yLabels[names(yLabels)=="Order"])))) %>%
      as.data.table()
    
    # Set up plot title
    pTitle <- paste0("Class: ", classRank, sep="")
    if(any(str_detect(tFS$rank, "Order"))){
      pTitle <- paste0(pTitle, ", nOrd=", ceiling(mean(tFS[(rank=="Order")]$nSum)), sep="")
    }
    if(any(str_detect(tFS$rank, "Family"))){
      pTitle <- paste0(pTitle, ", nFam=", ceiling(mean(tFS[(rank=="Family")]$nSum)), sep="")
    }
    if(any(str_detect(tFS$rank, "Genus"))){
      pTitle <- paste0(pTitle, ", nGen=", ceiling(mean(tFS[(rank=="Genus")]$nSum)), sep="")
    }
    if(any(str_detect(tFS$rank, "Species"))){
      pTitle <- paste0(pTitle, ", nSp=", ceiling(mean(tFS[(rank=="Species")]$nSum)), sep="")
    }
    
    # Create a ggplot bar plot
    ggtFail <- ggplot(tFS, aes(x = as.ordered(`Rank & Fold`), y = `tFail%`, 
                               fill = as.factor(klen))) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "Rank & Fold", 
           y = "Percentage of taxonomic groups where model training failed") +
      coord_flip() +
      theme_minimal() +
      theme(plot.title = element_text(size = 25, hjust=0.5),
            legend.position = 'top',
            legend.key.size = unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=18),
            axis.title = element_text(size = (22)),
            axis.text = element_text(size = (18))) +
      ggtitle(pTitle)
    
    # Write out to a file in the dirP directory
    ggsave(paste0(dirP, "/Barplot_tFail.png", sep=""), 
           plot=ggtFail,
           width = 12, height = 17, dpi = dpi, bg="white")
    
    # File names for bar plot
    plotNames[["TrainFail%"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "tFail")], sep="")
    
    # Rm unneeded vars
    rm(tFail, kSuccess, tFS, yLabels, ggtFail)
  }
  
  ### Bar plot of number of records (per fold/rank/klen) failing during classification ###
  
  # If groups were found to possess no discrim kmers
  if(nrow(rbindlist(modelResults[["ClassifyFail"]]))>0){
    # Unique folds present in classify fail
    uniqFolds <- sort(unique(rbindlist(modelResults[["ClassifyFail"]])$fold))
    
    # Count the total number of sequences per fold
    nSeqsFail <- rbindlist(modelResults[["ClassifyFail"]]) %>%
      select(c(recordID, fold)) %>%
      dplyr::distinct() %>%
      as.data.table()
    nSeqsSuccess <- rbindlist(modelResults[["CorrectIDs"]]) %>%
      select(c(recordID, fold)) %>%
      dplyr::distinct() %>%
      as.data.table()
    nSeqs <- rbind(nSeqsFail, nSeqsSuccess) %>%
      dplyr::distinct() %>%
      group_by(fold) %>%
      tally() %>%
      as.data.table()
    
    # Tally the number of records per rank/fold/klen where classifications failed
    cFail <- rbindlist(modelResults[["ClassifyFail"]]) %>% 
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nFail = n()) %>%
      mutate(`Rank & Fold` = paste0(rank, " ", fold, sep="")) %>%
      as.data.table()
    setkey(cFail, "rank","fold","klen")
    
    # Tally the number of records per rank/fold/klen where classifications succeeded
    cSuccess <- rbindlist(modelResults[["CorrectIDs"]], idcol = "rank") %>%
      select(c(rank, recordID, fold, klen)) %>%
      group_by(rank, klen, fold) %>%
      dplyr::summarise(nSuccess = n())  %>%
      as.data.table()
    setkey(cSuccess, "rank","fold","klen")
    
    # Merge cFail with cSuccess
    cFS <- merge(cFail, cSuccess) %>%
      rowwise() %>%
      mutate(`cFail%` = (nFail / (nFail + nSuccess)) * 100) %>%
      mutate(`Rank & Fold` = paste0(rank, " ", fold, sep=""))
    
    # Setting up y-axis label names for plot
    yLabels <- setNames(cFS$`Rank & Fold`, cFS$rank)
    
    # Set the labels as a factor for proper ordering by rank
    cFS <- cFS %>% 
      mutate(`Rank & Fold` = factor(`Rank & Fold`, 
                                    levels = c(unique(yLabels[names(yLabels)=="Species"]),
                                               unique(yLabels[names(yLabels)=="Genus"]), 
                                               unique(yLabels[names(yLabels)=="Family"]), 
                                               unique(yLabels[names(yLabels)=="Order"])))) %>%
      as.data.table()
    
    # Set up plot title
    pTitle <- paste0("Class: ", classRank, sep="")
    foreach(i=1:length(uniqFolds)) %do% {
      if(any(str_detect(cFS$fold, uniqFolds[i]))){
        pTitle <- paste0(pTitle, ", n", uniqFolds[i],"=", nSeqs[(fold==uniqFolds[i])]$n, sep="")
      }
    }
    pTitle <- pTitle[length(pTitle)]
    
    # Create a ggplot bar plot
    ggcFail <- ggplot(cFS, aes(x = as.ordered(`Rank & Fold`), y = `cFail%`, 
                               fill = as.factor(klen))) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "Rank & Fold", 
           y = "Percentage of test sequences where classifications failed") +
      coord_flip() +
      theme_minimal() +
      theme(plot.title = element_text(size = 25, hjust=0.5),
            legend.position = 'top',
            legend.key.size = unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=18),
            axis.title = element_text(size = (22)),
            axis.text = element_text(size = (18))) +
      ggtitle(pTitle) + 
      guides(fill = guide_legend(nrow = 1))
    
    # Write out to a file in the dirP directory
    ggsave(paste0(dirP, "/Barplot_cFail.png", sep=""), 
           plot=ggcFail,
           width = 12, height = 17, dpi = dpi, bg="white")
    
    # File names for bar plot
    plotNames[["ClassifyFail%"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "cFail")], sep="")
    
    # Rm unneeded vars
    rm(nSeqs, nSeqsFail, nSeqsSuccess, cFail, cSuccess, cFS, yLabels, ggcFail)
  }
  
  ### Bar plots of benchmarking results ###
  
  # If runBenchmark was selected as TRUE
  if(runBenchmark == TRUE){
    
    # Message update
    cat(green(paste0("Writing out benchmarking barplots to ", dirP, "...\n", sep="")))
    
    # Create naming vectors for the taxonomy column of each kData dataframe depending on numFolds & subclassRank
    ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
    if(numFolds >= 3){
      nameVecFolds <- paste0("Fold", 1:numFolds, sep="")
      nameVecTax <- c("recordID", ids, nameVecFolds)
    } else {
      nameVecTax <- c("recordID", ids)
    }
    
    # Create a dictionary with recordID and taxonomicIDs that can be referenced to determine
    # how many records belong to each fold
    recordFoldDict <- dfTaxa %>% 
      select(recordID, taxonomy) %>%
      separate(taxonomy, nameVecTax, ";")
    
    # Cols to be removed 
    rmCols <- c("ordID", "famID", "genID", "spID")
    rmCols <- append(rmCols, paste0("Fold", 1:numFolds, sep=""))
    
    # Grab all fold data for recordFoldDict and split by train, 
    # count number of records per fold
    recordTrain <- rbindlist(foreach(i=1:numFolds) %do% {
      # Subset the data into different folds
      recordFoldDict[which(recordFoldDict[,i+5] == "train"),] %>%
        mutate(folds = paste0("f", i))
    })
    # Take mean number of sequences across folds
    recordTrain <- recordTrain %>% 
      select(!all_of(rmCols)) %>% 
      group_by(folds) %>% 
      dplyr::summarise(n = n()) %>%
      dplyr::summarise(meanCnt = ceiling(mean(n)))
    
    # Remove unneeded vars
    rm(recordFoldDict, rmCols, ids)
    
    # Read in bench data
    kfind <- read_dbTable(dbPath, "bench_kmer_Find")
    
    # Sum up data from kFind function
    kfindSum <- kfind %>%
      group_by(klen) %>%
      mutate(sumTime = sum(`elapsedTime (mins)`)) %>%
      mutate(sumMem = sum(`memUsage (Mb)`)) %>%
      select(c(klen, sumTime, sumMem)) %>%
      dplyr::distinct() %>%
      arrange(klen)
    
    # Modify kFind with kFindSum values
    kfind <- kfind %>% 
      filter(subFuncName == "kPos") %>%
      select(!subFuncName) %>% 
      mutate(`elapsedTime (mins)` = kfindSum$sumTime) %>%
      mutate(`memUsage (Mb)` = kfindSum$sumMem)
    
    # Remove unneeded vars
    rm(kfindSum)
    
    # Read in the benchmarking data for the rest of the functions
    kdiscrim <- read_dbTable(dbPath, "bench_kmer_discrimFind") %>% select(!subFuncName) 
    ktrain <- read_dbTable(dbPath, "bench_kphmm_Train") %>% select(!subFuncName)
    ktest <- read_dbTable(dbPath, "bench_kphmm_ClassifyCv")
    
    # Divide ktest into per sequence benchmarks and total benchmarks
    
    # Total
    ktestTotal <- ktest %>%
      mutate(`elapsedTime (mins)` = `elapsedTime Total (mins)`) %>%
      mutate(`memUsage (Mb)` = `memUsage Total (Mb)`) %>%
      mutate(`dbUsage (Mb)` = 0) %>%
      select(!c(`elapsedTime Total (mins)`, 
                `memUsage Total (Mb)`,
                `elapsedTime AvgPerSeq (mins)`,
                `memUsage AvgPerSeq (Mb)`,
                nSeqs))
    
    # Per sequence
    ktestInd <- ktest %>%
      mutate(`elapsedTime (mins)` = `elapsedTime AvgPerSeq (mins)`) %>%
      mutate(`memUsage (Mb)` = `memUsage AvgPerSeq (Mb)`) %>%
      mutate(`dbUsage (Mb)` = 0) %>%
      select(!c(`elapsedTime Total (mins)`, 
                `memUsage Total (Mb)`,
                `elapsedTime AvgPerSeq (mins)`,
                `memUsage AvgPerSeq (Mb)`,
                nSeqs))
    
    # Find per sequence benchmarks for ktrain
    ktrainInd <- ktrain %>%
      mutate(nSeqsUsed = recordTrain$meanCnt[1]) %>%
      rowwise() %>%
      mutate(`memUsage (Mb)` = `memUsage (Mb)` / nSeqsUsed) %>%
      mutate(`elapsedTime (mins)` = `elapsedTime (mins)` / nSeqsUsed)
    
    # Total elapsedTime and memory used for all functions
    benchTotal <- rbind(kfind, kdiscrim, ktrain, ktestTotal) %>%
      mutate(funcName = factor(funcName, 
                               levels = c("kmer_Find",
                                          "kmer_discrimFind", 
                                          "kphmm_Train", 
                                          "kphmm_ClassifyCV")))
    
    # elapsed time and memory usage per sequence for kphmm_Train and kphmm_ClassifyCV
    benchInd <- rbind(ktrainInd, ktestInd)
    
    # Plot for total time taken
    gBenchTotalTime <- ggplot(benchTotal, aes(x = as.factor(klen), y = `elapsedTime (mins)`, fill = klen)) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "", y = "Average elapsed time (mins)") +
      theme(plot.title = element_text(size=25, hjust=0.5),
            legend.position='top',
            legend.key.size=unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=18),
            axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.x = element_text(size = 22)) + 
      guides(fill = guide_legend(nrow = 1)) + 
      facet_wrap(~ funcName, scales = "free_x", nrow = 4, ncol = 1) +
      coord_flip()
    
    # Plot for total memory usage
    gBenchTotalMem <- ggplot(benchTotal, aes(x = as.factor(klen), y = `memUsage (Mb)`, fill = klen)) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "", y = "Average memory usage (Mb)") +
      theme(plot.title = element_text(size=25, hjust=0.5),
            legend.position='top',
            legend.key.size=unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=18),
            axis.title=element_text(size=22),
            axis.text=element_text(size=22),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.x = element_text(size = 22)) + 
      guides(fill = guide_legend(nrow = 1)) + 
      facet_wrap(~ funcName, scales = "free_x", nrow = 4, ncol = 1) +
      coord_flip()
    
    # Use ggarrange to combine into one plot
    gBenchTotal <- ggarrange(gBenchTotalTime, gBenchTotalMem, ncol=2, nrow=1, common.legend = TRUE, legend="top") 
    gBenchTotal <- annotate_figure(gBenchTotal, top = text_grob(paste0("Class: ", classRank, ", nFolds=", numFolds, 
                                                                       ", nCores=", numCores, ", nTest=", ktest$nSeqs[1], 
                                                                       ", nTrain=",  recordTrain$meanCnt[1], sep=""), size = 25))
    # Save out to a file
    ggsave(paste0(dirP, "/Barplot_BenchTotals.png", sep=""), 
           plot=gBenchTotal,
           width = 12, height = 20, dpi = dpi, bg="white")
    
    # Plot for total time taken
    gBenchIndTime <- ggplot(benchInd, aes(x = as.factor(klen), y = `elapsedTime (mins)`, fill = klen)) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "", y = "Average elapsed time per sequence (mins)") +
      theme(plot.title = element_text(size=25, hjust=0.5),
            legend.position='top',
            legend.key.size=unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=19),
            axis.title=element_text(size=19),
            axis.text=element_text(size=19),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.x = element_text(size = 22)) + 
      guides(fill = guide_legend(nrow = 1)) + 
      facet_wrap(~ funcName, scales = "free_x", nrow = 2, ncol = 1) +
      coord_flip()
    
    # Plot for total memory usage
    gBenchIndMem <- ggplot(benchInd, aes(x = as.factor(klen), y = `memUsage (Mb)`, fill = klen)) +
      geom_bar(stat = "Identity", position = position_dodge2(preserve = "single"), colour = "black") +
      scale_fill_viridis(discrete = TRUE, name = "kLength") +
      labs(x = "", y = "Average memory usage per sequence (Mb)") +
      theme(plot.title = element_text(size=25, hjust=0.5),
            legend.position='top',
            legend.key.size=unit(1, 'cm'),
            legend.title=element_text(size=22),
            legend.text=element_text(size=19),
            axis.title=element_text(size=19),
            axis.text=element_text(size=19),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.x = element_text(size = 22)) + 
      guides(fill = guide_legend(nrow = 1)) + 
      facet_wrap(~ funcName, scales = "free_x", nrow = 2, ncol = 1) +
      coord_flip()
    
    # Use ggarrange to combine into one plot
    gBenchInd <- ggarrange(gBenchIndTime, gBenchIndMem, ncol=2, nrow=1, common.legend = TRUE, legend="top")
    gBenchInd <- annotate_figure(gBenchInd, top = text_grob(paste0("Class: ", classRank, ", nFolds=", numFolds, 
                                                                   ", nCores=", numCores, ", nTest=", ktest$nSeqs[1], 
                                                                   ", nTrain=",  recordTrain$meanCnt[1], sep=""), size = 25))
    # Save out to a file
    ggsave(paste0(dirP, "/Barplot_BenchPerSeq.png", sep=""), 
           plot=gBenchInd,
           width = 12, height = 15, dpi = dpi, bg="white")
    
    # File names for bar plot
    plotNames[["BenchTotals"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "BenchTotals")], sep="")
    plotNames[["BenchperSeq"]] <- paste0(dirP, "/", list.files(dirP)[str_detect(list.files(dirP), "BenchPerSeq")], sep="")
    
    # Rm unneeded vars
    suppressWarnings(rm(gBenchIndMem, gBenchIndTime, gBenchInd, gBenchTotal, 
                        gBenchTotalTime, gBenchTotalMem, ktestTotal, kfind, benchInd, benchTotal,
                        kdiscrim, ktrain, ktrainInd, ktest, ktestInd, ktestTotal))
  }
  
  # Return the plot name from the function
  return(plotNames)
}