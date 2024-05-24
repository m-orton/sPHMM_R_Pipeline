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
  classPos <- lapply(dnaNames, function(x){ which(x == classRank) })
  suppressWarnings(taxaData <- foreach(i=1:length(dnaNames)) %do% dnaNames[[i]][classPos[[i]]:(classPos[[i]]+5)])
  
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
  
  # Success message
  cat(green(paste0("Input fasta file successfully loaded for class ", classRank, ".\n")))
  
  # remove all unecessary vars created
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
seq_Filter <- function(dfTaxa, nFilter, gFilter){
  dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "N|n") / nchar(dfTaxa$dna)) <= nFilter),]
  dfTaxa <- dfTaxa[which((str_count(dfTaxa$dna, "-") / nchar(dfTaxa$dna)) <= gFilter),]
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
  dir.create(dir, showWarnings = FALSE)
  
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
      folds <- unite(folds, col='folds', colnames(folds), sep=';')
      
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
      folds <- unite(folds, col='folds', colnames(folds), sep=';')
      
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

##### Discriminitive Kmer Functions #####

# Internal function to load the C++ kmer counting script for grabbing the positional data of the kmers (using the Rcpp R package)
sourceCpp("kmercntwritecolumns.cpp")

# The main function used for kmer extraction function
kmer_Find <- function(dfTaxa, dir, klen){
  
  # Creating directory for kmer files
  dirK <- paste0(dir, "/KmerData")
  dir.create(dirK, showWarnings = FALSE)
  
  # Convert to DNAbin format
  dnaBin <- as.DNAbin(DNAStringSet(dfTaxa$dna))
  names(dnaBin) <- dfTaxa$taxonomy
  
  # Extract kmers looping through each klen
  kmerFiles_c <- replicate(length(klen), list())
  kmerFiles_p <- replicate(length(klen), list())
  
  foreach(i=1:length(klen)) %do% {
    # Update progress  
    cat(green(paste0("Extracting kmers", " (klen=", klen[i], ")...\n", sep="")))
    
    # Create filenames for count matrices
    kmerFiles_c[[i]] <- paste(dirK, "/raw_klen_", klen[i], ".tsv", sep="")
    
    # Create filenames for position data
    kmerFiles_p[[i]] <- paste(dirK, "/raw_klen_pos", klen[i], ".tsv", sep="")
    
    # use kmer package to grab kmer data for all sequences
    count_matrix <- kmer::kcount(dnaBin, k = klen[i], residues = "DNA")
    
    # Sum the counts across columns for with Rfast package
    col_totals <- Rfast::colsums(count_matrix)
    
    # Calculate frequencies for columns
    count_matrix <- count_matrix / col_totals
    
    # Convert to dataframe format
    count_matrix <- as.data.frame(count_matrix) %>% 
      rownames_to_column(var="taxonomy")
    
    # Write count data to tsv format
    write_tsv(count_matrix, kmerFiles_c[[i]])
    
    # Write position data to tsv format with c++ function for each klen
    kmer_countToFile(as.character(dfTaxa$dna), dfTaxa$taxonomy, klen[i], kmerFiles_p[[i]])
    
    # remove count matrix and totals each iteration after writing to disk
    rm(count_matrix, col_totals)
  }
  
  # Remove dnaBin
  rm(dnaBin)
  
  # Name files
  names(kmerFiles_c) <- paste0("klen", klen, sep="")
  names(kmerFiles_p) <- paste0("klen", klen, sep="")
  kmerFiles <- list("CountData"=kmerFiles_c, "PosData"=kmerFiles_p)

  # Final message
  cat(green(paste0("Kmers (klen=", paste(klen, collapse=','), ") have been extracted and written to the ", dir, 
                   "/kmerData folder in tsv format.\n")))
  
  # return files created for each kmer length
  return(kmerFiles)
}

# New function for finding discriminative kmers using the doFuture package for parallel processing
kmer_discrimFind <- function(kmerFiles, dfTaxa, numFolds, subclassRank, klen, sigVal, numMinK, numMaxK, runParallel, numCores){
  
  # List for kData
  kData <- replicate(length(klen), list())
  
  # List for wilcox test results data
  wTest <- replicate(length(klen), list())
  
  # Create naming vectors for the taxonomy column of each kData dataframe depending on numFolds & subclassRank
  ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
  if(numFolds >= 3){
    nameVecFolds <- paste0("Fold", 1:numFolds, sep="")
    nameVecTax <- c("recordID", ids, nameVecFolds)
  } else {
    nameVecTax <- c("recordID", ids)
  }
  
  # Iterate through each kmer length extracting the results
  foreach(i=1:length(klen)) %do% {
    # Read in raw kmer data
    kData[[i]] <- read_tsv(as.character(kmerFiles[[1]][[i]]), show_col_types = FALSE) %>%
      separate(taxonomy, nameVecTax, sep = ";") 
    
    # Perform wilcox rank sum tests depending on which taxonomic ranks are being run
    # Order level
    if("Order" %in% subclassRank){
      rank <- "Order"
      message('\r', green(paste0("Searching kmers... klen=", klen[i], ", Rank: Order     ", sep="")), appendLF = FALSE)
      
      # Find all unique orders
      uniqOrd <- unique(kData[[i]]$ordID)
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # If parallel processing is used
      if(runParallel==TRUE){
        #registerDoFuture()
        plan(multisession, workers=numCores)
        # Iterate through each unique order
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqOrd)) %dofuture% {
          .kmer_wilcoxTest(kData[[i]], uniqOrd[j], rank, numFolds, sigVal)
        }))
        plan(sequential)
        # If no parallel processing is being performed
      } else {
        # Iterate through each unique order
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqOrd)) %do% {
          .kmer_wilcoxTest(kData[[i]], uniqOrd[j], rank, numFolds, sigVal)
        }))
      }
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK)
    }
    # Family Level
    if("Family" %in% subclassRank){
      rank <- "Family"
      message('\r', green(paste0("Searching kmers... klen=", klen[i], ", Rank: Family     ", sep="")), appendLF = FALSE)
      
      # Find all unique orders
      uniqFam <- unique(kData[[i]]$famID)
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # If parallel processing is used
      if(runParallel==TRUE){
        #registerDoFuture()
        plan(multisession, workers=numCores)
        # Iterate through each unique family
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqFam)) %dofuture% {
          .kmer_wilcoxTest(kData[[i]], uniqFam[j], rank, numFolds, sigVal)
        }))
        plan(sequential)
        # If no parallel processing is being performed
      } else {
        # Iterate through each unique family
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqFam)) %do% {
          .kmer_wilcoxTest(kData[[i]], uniqFam[j], rank, numFolds, sigVal)
        }))
      }
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK)
    }
    # Genus Level
    if("Genus" %in% subclassRank){
      rank <- "Genus"
      message('\r', green(paste0("Searching kmers... klen=", klen[i], ", Rank: Genus    ", sep="")), appendLF = FALSE)
      
      # Find all unique orders
      uniqGen <- unique(kData[[i]]$genID)
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # If parallel processing is used
      if(runParallel==TRUE){
        #registerDoFuture()
        plan(multisession, workers=numCores)
        # Iterate through each unique genus
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqGen)) %dofuture% {
          .kmer_wilcoxTest(kData[[i]], uniqGen[j], rank, numFolds, sigVal)
        }))
        plan(sequential)
        # If no parallel processing is being performed
      } else {
        # Iterate through each unique genus
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqGen)) %do% {
          .kmer_wilcoxTest(kData[[i]], uniqGen[j], rank, numFolds, sigVal)
        }))
      }
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK)
    }
    # Species Level
    if("Species" %in% subclassRank){
      rank <- "Species"
      message('\r', green(paste0("Searching kmers... klen=", klen[i], ", Rank: Species    ", sep="")), appendLF = FALSE)
      
      # Find all unique orders
      uniqSp <- unique(kData[[i]]$spID)
      
      # Sublist for wTest list
      wTest[[i]][[rank]] <- list()
      
      # If parallel processing is used
      if(runParallel==TRUE){
        #registerDoFuture()
        plan(multisession, workers=numCores)
        # Iterate through each unique species
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqSp)) %dofuture% {
          .kmer_wilcoxTest(kData[[i]], uniqSp[j], rank, numFolds, sigVal)
        }))
        plan(sequential)
        # If no parallel processing is being performed
      } else {
        # Iterate through each unique species
        wTest[[i]][[rank]] <- suppressWarnings(rbindlist(foreach(j=1:length(uniqSp)) %do% {
          .kmer_wilcoxTest(kData[[i]], uniqSp[j], rank, numFolds, sigVal)
        }))
      }
      # Filter by numMinK and numMaxK
      wTest[[i]][[rank]] <- .kmer_Filter(wTest[[i]][[rank]], numFolds, numMinK, numMaxK)
    }
  }
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
  
  # Remove columns not needed
  wTest <- wTest %>% select(!c("kCount_Group"))
  
  # kmerResults list to return from the function
  kmerResults <- list("kResults" = wTest,
                      "kTally" = tallyK)
  
  # Remove unneeded vars
  rm(wTest, tallyK, kData)
  
  # Return kmerResults
  return(kmerResults)
}

# Internal function to perform wilcoxon ranked sum tests

# For each kmer group and taxonomic group at a specified taxonomic rank, perform a series of 
# group vs non-group comparisons. In this case group refers to kmer frequencies of target taxonomic group
# and non-group is any other kmer frequencies not belonging to that group (ex: all Bovidae kmer freqs vs all Non-Bovidae kmer freqs)
# A two-sided, non-paired wilcoxon rank-sum test is then performed for each group/nongroup and p-values are extracted from the test
.kmer_wilcoxTest <- function(kData, group, rank, numFolds, sigVal){
  
  # Subset data depending on which rank is being run
  if(rank == "Order"){
    col_g <- kData %>% filter(ordID == group)
    col_ng <- kData %>% filter(ordID != group)
  } else if (rank == "Family"){
    col_g <- kData %>% filter(famID == group)
    col_ng <- kData %>% filter(famID != group)
  } else if (rank == "Genus"){
    col_g <- kData %>% filter(genID == group)
    col_ng <- kData %>% filter(genID != group)
  } else {
    col_g <- kData %>% filter(spID == group)
    col_ng <- kData %>% filter(spID != group)
  }
  
  # If cross-val tests are being performed, create the necessary fold data first before proceeding
  if(numFolds < 3){
    if(!is.null(ncol(col_g))){
      wilcox <- col_wilcoxon_twosample(
        as.matrix(col_g[,6:ncol(col_g)]),
        as.matrix(col_ng[,6:ncol(col_ng)]),
        null = 0,
        alternative = "two.sided",
        exact = NA,
        correct = TRUE) %>% 
        rownames_to_column(var = "kmers") %>%
        adjust_pvalue(method = "bonferroni") %>%
        mutate(group = group) %>% 
        mutate(rank = rank) %>% 
        select(group, rank, kmers, pvalue.adj) %>%
        filter(pvalue.adj < sigVal & pvalue.adj > 0)
    } else {
    }
  } else {
    positions <- c(1,1:11)
    wilcox <- replicate(numFolds, list())
    # Find the rows labelled as train from each fold
    foreach(i=1:numFolds) %do% {
      if(!is.null(ncol(col_g))){
        wilcox[[i]] <- col_wilcoxon_twosample(
          as.matrix(col_g %>% filter(select(.,i+5) == "test") %>% select(!positions)),
          as.matrix(col_ng %>% filter(select(.,i+5) == "test") %>% select(!positions)),
          null = 0,
          alternative = "two.sided",
          exact = NA,
          correct = TRUE) %>% 
          rownames_to_column(var = "kmers") %>%
          adjust_pvalue(method = "bonferroni") %>%
          mutate(group = group) %>% 
          mutate(rank = rank) %>% 
          mutate(fold = paste("f", i, sep="")) %>%
          select(group, rank, kmers, fold, pvalue.adj) %>%
          filter(pvalue.adj < sigVal & pvalue.adj > 0)
      } else {
      }
    }
    # rbind wilcox test results across folds
    wilcox <- rbindlist(wilcox)
  }
  
  # remove unneeded vars
  rm(col_g, col_ng, positions)
  
  # return test results
  return(wilcox)
}

# Internal function to filter discriminative kmers by numMaxK and numMinK
.kmer_Filter <- function(wTest, numFolds, numMinK, numMaxK){
  
  # Number of unique groups present
  uniqGroups <- length(unique(wTest$group))
  
  if(numFolds >=3){
    # Determine kCoverage_Intergroup (proportion of all unique groups that possess said kmer) and filter by kCount_Group >= numMinK
    wTest <- wTest %>% 
      group_by(fold, kmers) %>% 
      add_count(kmers, name = "kCoverage_Intergroup") %>%
      mutate(kCoverage_Intergroup = kCoverage_Intergroup / uniqGroups) %>%
      ungroup() %>%
      group_by(fold, group) %>%
      add_count(group, name = "kCount_Group") %>%
      filter(kCount_Group >= numMinK)
    
    # Arrange by both by numMaxK if value is not NA
    if(!is.na(numMaxK)){
    wTest <- wTest %>%
      # Attempt to sort by lowest p-value & kCoverage_Intergroup 
      arrange(kCoverage_Intergroup, pvalue.adj) %>%
      slice_head(n = numMaxK)
    }
  } else {
    # Determine kCoverage_Intergroup (proportion of all unique groups that possess said kmer) and filter by kCount_Group >= numMinK
    wTest <- wTest %>% 
      group_by(kmers) %>% 
      add_count(kmers, name = "kCoverage_Intergroup") %>%
      mutate(kCoverage_Intergroup = kCoverage_Intergroup / uniqGroups) %>%
      ungroup() %>%
      group_by(group) %>%
      add_count(group, name = "kCount_Group") %>%
      filter(kCount_Group >= numMinK)
    
    # Arrange by both by numMaxK if value is not NA
    if(!is.na(numMaxK)){
      wTest <- wTest %>%
        # Attempt to sort by lowest p-value & kCoverage_Intergroup 
        arrange(kCoverage_Intergroup, pvalue.adj) %>%
        slice_head(n = numMaxK)
    }
  }
  # Return wTest
  return(wTest)
}

##### PHHM Model Setup and Training Functions #####

# Internal function to actually construct the kPHMM itself when using fold data for cross-validation
.kphmm_runTrainF <- function(kMat, weights, residue, rankName, groupName, klenList, foldName, numMaxiter, deltaLL, numCores) {
  
  # Messages for success or failure with tryCatch
  success <- paste0("Success! kPHHM trained for group: ", groupName, "                           ", sep="")
  warning <- paste0("Warning: The kPHHM for group: ", groupName, " could not be constructed since the function kphmm_runTrainF failed.", 
                    "                           ", sep="")
  
  # Try catch handling for running the PHMM itself with the function derive.PHMM
  tryCatch(
    # If the derivePHHM model is successful
    { 
      construct.PHMM <- derivePHMM(kMat, seqweights = as.numeric(weights), residues = residue, name = groupName, refine = "Baumwelch", 
                                   pseudocounts = "background", maxiter = numMaxiter, deltaLL = deltaLL, cores = numCores, quiet = TRUE)
      construct.PHMM <- list("Group"= groupName, "Rank"= rankName, "klen"=klenList, "fold"=foldName, "model"=construct.PHMM, "messages"=success)
      # success message
      message('\r', green(success), appendLF = FALSE)
      Sys.sleep(0.01)
      return(construct.PHMM)
    },
    # If the derivePHHM model is not successful
    error = function(e) {
      construct.PHMM <- list("Group"=groupName, "Rank"= rankName, "klen"=klenList, "fold"=foldName, "model"=NULL,"messages"=warning)
      # failure message
      message('\r', red(warning), appendLF = FALSE)
      Sys.sleep(0.01)
      return(construct.PHMM)
    }
  )
}

# Internal function used to actually construct the kPHMM itself when not using fold data for cross-validation
.kphmm_runTrainNF <- function(kMat, weights, residue, rankName, groupName, klenList, numMaxiter, deltaLL, numCores) {
  
  # Messages for success or failure with tryCatch
  success <- paste0("Success! kPHHM trained for group: ", groupName, "                           ", sep="")
  warning <- paste0("Warning: The kPHHM for group: ", groupName, " could not be constructed since the function kphmm_runTrainNF failed.",
                    "                           ", sep="")
  
  # Try catch handling for running the PHMM itself with the function derive.PHMM
  tryCatch(
    # If the derivePHHM model is successful
    { 
      construct.PHMM <- derivePHMM(kMat, seqweights = as.numeric(weights), residues = residue, name = groupName, refine = "Baumwelch", 
                                   pseudocounts = "background", maxiter = numMaxiter, deltaLL = deltaLL, cores = numCores, quiet = TRUE)
      construct.PHMM <- list("Group"= groupName, "Rank"= rankName, "klen"=klenList, "model"=construct.PHMM, "messages"=success)
      # success message
      message('\r', green(success), appendLF = FALSE)
      Sys.sleep(0.01)
      return(construct.PHMM)
    },
    # If the derivePHHM model is not successful
    error = function(e) {
      construct.PHMM <- list("Group"=groupName, "Rank"= rankName, "klen"=klenList, "model"=NULL,"messages"=warning)
      # failure message
      message('\r', red(warning), appendLF = FALSE)
      Sys.sleep(0.01)
      return(construct.PHMM)
    }
  )
}

# Main function used to do all setup & training of kPHMM models (when performing cross-validation) using the test labeled sequences from each fold or
# if not performing cross-validation tests, the training sequences provided by the user
kphmm_Train <- function(kmerResults, kmerFiles, klen, numFolds, chunkSize, numMaxKPos, numMinKSeq, gapThreshold, numMaxiter, subclassRank, deltaLL, numCores){
  
  # Starting message
  cat(green(paste0("Organizing kmer data for model training (this may take a few moments)...\n")))
  
  # Vector for subclassRank column ids
  ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
  
  if(numFolds >= 3){
    # Naming vectors for taxonomy and fold data
    nameVecFolds <- paste0("f", 1:numFolds, sep="")
    nameVecTax <- c("recordID", ids, nameVecFolds)
    
    # Grab kmer data from kmerResults
    kmers <- kmerResults[['kResults']] %>% 
      select(fold, kmers, group, klen, rank)
    
    # List for kPosData
    kPosData <- list()
    
    # Grab the position data from the raw kmer data files
    foreach(i=1:length(klen)) %do% {
      # Read the data in chunkwise
      kPosData[[i]] <- read_chunkwise(as.character(kmerFiles[[2]][[i]]), sep="\t", chunk_size = 100000)
      # select position, kmers and taxonomy columns and collect from chunkwise objects
      kPosData[[i]] <- kPosData[[i]] %>%
        select(cols = c(3, 4, 6)) %>%
        # rename pos, kmers, taxonomy columns
        mutate(pos = cols1) %>%
        mutate(kmers = cols2) %>%
        mutate(taxonomy = cols3) %>%
        # select all relevant columns and collect data
        select(pos, kmers, taxonomy) %>%
        collect() %>% 
        separate(taxonomy, nameVecTax, sep = ";") %>%
        mutate(klen = klen[i])
    }
    
    # Rbind positional data across klengths
    kPosData <- data.frame(rbindlist(kPosData))
    
    # list for folds
    folds <- replicate(numFolds, list())
    
    # Find the rows labelled as train from each fold
    foreach(i=1:numFolds) %do% {
      # Subset the data into different folds
      folds[[i]] <- kPosData[which(kPosData[,i+7] == "train"),]
      # modify kmers column to also include fold number
      folds[[i]]$folds <- paste0("f", i)
    }
    
    # Rbind all folds together into one large dataframe and remove duplicates
    kPosData <- rbindlist(folds) %>% 
      select(folds, pos, kmers, recordID, ordID, famID, genID, spID, klen) %>%
      distinct()
    
    # Divide the dataframe by subclassRank and reorganize columns for rank and group
    # Order
    kPosData_o <- kPosData %>% 
      select(folds, pos, kmers, recordID, ordID, klen) %>% 
      mutate(group = ordID) %>% 
      mutate(rank = "Order") %>% 
      select(!ordID)
    
    # Family
    kPosData_f <- kPosData %>% 
      select(folds, pos, kmers, recordID, famID, klen) %>% 
      mutate(group = famID) %>% 
      mutate(rank = "Family") %>% 
      select(!famID)
    
    # Genus
    kPosData_g <- kPosData %>% 
      select(folds, pos, kmers, recordID, genID, klen) %>% 
      mutate(group = genID) %>% 
      mutate(rank = "Genus") %>% 
      select(!genID)
    
    # Species
    kPosData_s <- kPosData %>% 
      select(folds, pos, kmers, recordID, spID, klen) %>% 
      mutate(group = spID) %>% 
      mutate(rank = "Species") %>% 
      select(!spID)
    
    # Combine back together with rbind
    kPosData <- rbind(kPosData_o, kPosData_f, kPosData_g, kPosData_s)
    
    # Remove unneeded vars
    rm(folds, kPosData_o, kPosData_f, kPosData_g, kPosData_s)
    
    # Ensure both kmers and kPosData are in dataframe format
    kmers <- data.frame(kmers)
    kPosData <- data.frame(kPosData)
    
    # Compare groups here with groups and folds in the kmers dataframe and subset out those that are not present (since those are unusable)
    kmers$group_fold <- paste0(kmers$folds, "_", kmers$group, "_", kmers$klen, sep="")
    kPosData$group_fold <- paste0(kPosData$folds, "_", kPosData$group, "_", kPosData$klen, sep="")
    kPosData <- kPosData[kPosData$group_fold %in% kmers$group_fold,]
    
    # Split into a list by folds, group & klen for both kPosData & kmers dataframes
    
    # kPosData dataframe
    kPosData <- kPosData %>% 
      select(!group_fold) %>% 
      group_by(folds, group, klen) %>%
      group_split() 
    # kPosData <- kPosData %>% 
    #  group_by(group_fold) %>%
    #  group_split()
    
    # kmers dataframe
    kmers <- kmers %>%
     select(!group_fold) %>%
     group_by(folds, group, klen) %>%
     group_split()
    # kmers <- kmers %>%
    #  group_by(group_fold) %>%
    #  group_split()
    
    # This set of for loops will convert to DT for every group
    kmers <- foreach(i=1:length(kmers)) %do% data.table(kmers[[i]])
    kPosData <- foreach(i=1:length(kPosData)) %do% data.table(kPosData[[i]])
    
    # Setkeys for every dataframe across kPosData and kmers
    foreach(i=1:length(kmers)) %do% {
      setkey(kPosData[[i]], "kmers")
      setkey(kmers[[i]], "kmers")
    }
    
    # Merge kmer dataframes with kPosData dataframes using kmers as the key
    kPosData <- foreach(i=1:length(kPosData)) %do% merge.data.table(kPosData[[i]], kmers[[i]])
    
    # Find all unique kmers 
    uniqKmers <- foreach(i=1:length(kmers)) %do% unique(kmers[[i]]$kmers)
    
    # Count number of recordIDs
    recordNum <- foreach(i=1:length(kPosData)) %do% length(unique(kPosData[[i]]$recordID))
    
    # Tally up the number of disrim kmers at each position / total number of sequences
    kTally <- list()
    foreach(i=1:length(kPosData)) %do% {
      kTally[[i]] <- kPosData[[i]] %>%
        group_by(pos, kmers) %>%
        tally() %>% 
        mutate(n = n / recordNum[[i]]) 
    }
    
    # Search for positions with the highest numbers of discrim kmers, pick the top number of positions (as defined by numMaxKPos)
    foreach(i=1:length(kTally)) %do% {
      kTally[[i]] <- subset(kTally[[i]], kTally[[i]]$kmers %in% uniqKmers[[i]]) %>%
        arrange(desc(n)) %>%
        ungroup(pos, kmers) %>%
        slice_head(n = numMaxKPos)
    }
    
    # Subset by those positions (highest numbers of discrim kmers), 
    kMat <- list()
    foreach(i=1:length(kTally)) %do% {
      kMat[[i]] <- subset(kPosData[[i]], kPosData[[i]]$pos %in% kTally[[i]]$pos) %>%
        mutate(data = paste0(rank.x, ";", group.x, ";", klen.x, ";", folds.x,  sep="")) %>%
        select(data, recordID, pos, kmers)
    }
    
    # Delete any duplicate rows (across all columns), relabel 'NC' record ids to be unique by rownumber and spread out data in a matrix like format for each list element, 
    # fill by gap symbol
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- dplyr::distinct(kMat[[i]]) %>%
        mutate(recordID = ifelse(recordID == "NC", paste0("NC", row_number(), sep=""), recordID)) %>%
        spread(pos, kmers, fill = "-")
    }
    
    # Find gappy kmer positions and filter them out, removal of these kmer positions is governed by the gapThreshold parameter
    kGap <- list()
    foreach(i=1:length(kMat)) %do% {
      kGap[[i]] <- kMat[[i]] %>%
        gather(pos, kmers) %>%
        filter(str_detect(pos,"^\\s*[0-9]*\\s*$")) %>%
        group_by(pos) %>% 
        dplyr::summarise(sumG = sum(kmers == "-") / nrow(kMat[[i]])) %>%
        filter(sumG <= gapThreshold)
    }
    
    # kmerCheck var to find which df's are empty from kGap
    kmerCheck <- foreach(i=1:length(kGap)) %do% {
      nrow(kGap[[i]]) == 0
    }
    
    # Subset kMat, kGap and uniqKmers by kmerCheck
    kMat <- kMat[!unlist(kmerCheck)]
    kGap <- kGap[!unlist(kmerCheck)]
    uniqKmers <- uniqKmers[!unlist(kmerCheck)]
    
    # Remove kmerCheck
    rm(kmerCheck)
    
    # Subset kMat by kGap positions
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- cbind.data.frame(kMat[[i]][,1:2], subset(kMat[[i]], select = colnames(kMat[[i]]) %in% kGap[[i]]$pos))
    }
    
    # Count how many discriminitive kmers each sequence has (identified by recordID) and filter by a minimum number of discrim kmers per recordID (as defined by numMinKSeq)
    kMat2 <- list()
    foreach(i=1:length(kMat)) %do% {
      kMat2[[i]] <- kMat[[i]] %>%
        mutate(count = across(.cols = 3:ncol(kMat[[i]]), ~ .x %in% uniqKmers[[i]])) %>%
        rowwise() %>%
        mutate(countDiscrim = across(.cols = contains("count"), .fns = sum)) %>%
        select(recordID, countDiscrim) %>%
        filter(countDiscrim > numMinKSeq[1])
    }
    
    # Add a weighted count column to assign weights to each sequence for model setup based on number of discrim kmer positions they posses
    # ex: a sequence with the all possible positions would be assigned the highest possible weight proportional to the total number of sequences
    # note: average weight across all sequences must equal 1 and sum of all weights must equal total number of sequences
    foreach(i=1:length(kMat2)) %do% {
      kMat2[[i]]$sumCount <- sum(kMat2[[i]]$countDiscrim)
      kMat2[[i]]$weight <- (kMat2[[i]]$countDiscrim / kMat2[[i]]$sumCount) * nrow(kMat2[[i]])
    }
    
    # Final subset by recordID to get the matrices that will be used in the PHHM's and add the previously calculated weight column
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- subset(kMat[[i]], kMat[[i]]$recordID %in% kMat2[[i]]$recordID)
      kMat[[i]]$weight <- kMat2[[i]]$weight
    }
    
    # Remove unneeded vars
    rm(kMat2, kGap, kTally, uniqKmers, recordNum, kPosData, kmers)
    
    # Transfer recordID column to rowname
    kMat <- foreach(i=1:length(kMat)) %do% {column_to_rownames(as.tibble(kMat[[i]]), var = "recordID")}
    
    # Convert kMat to matrix format
    foreach(i=1:length(kMat)) %do% {kMat[[i]] <- as.matrix(kMat[[i]])}
    
    # Check for groups that possess no discrim kmer positions (after the above steps) and filter them out
    kmerCheck <- foreach(i=1:length(kMat)) %do% {
      length(kMat[[i]]) == 0 | ncol(kMat[[i]]) < 3
    }
    
    # Subset kMat by kmerCheck
    kMat <- kMat[!unlist(kmerCheck)]
    
    # Remove unneeded vars
    rm(kmerCheck)
    
    # Take the first rowname element from each list element
    nameList <- unlist(foreach(i=1:length(kMat)) %do% {str_split(kMat[[i]][1,1], ";")}, recursive = FALSE)
    
    # rankName from nameList
    if(length(klen) != 1){
      rankName <- sapply(nameList, `[`, 1)
    } else {
      rankName <- as.list(sapply(nameList, `[`, 1))
    }
    
    # groupName from nameList
    if(length(klen) != 1){
      groupName <- sapply(nameList, `[`, 2)
    } else {
      groupName <- as.list(sapply(nameList, `[`, 2))
    }
    
    # klen from nameList
    if(length(klen) != 1){
      klenList <- sapply(nameList, `[`, 3)
    } else {
      klenList <- as.list(sapply(nameList, `[`, 3))
    }
    
    # foldName from nameList
    if(length(klen) != 1){
      foldName <- sapply(nameList, `[`, 4)
    } else {
      foldName <- as.list(sapply(nameList, `[`, 4))
    }
    
    # Extract the weight column
    weights <- foreach(j=1:length(kMat)) %do% kMat[[j]][, ncol(kMat[[j]])]
    
    # Remove the data column and the weight column from each matrix 
    kMat <- foreach(i=1:length(kMat)) %do% kMat[[i]][,-1]
    kMat <- foreach(i=1:length(kMat)) %do% kMat[[i]][,-ncol(kMat[[i]])]
    
    # Create vectors of residues for the PHMM for each klen (all possible kmers for a given klen using gtools::permutations)
    residues <- list()
    bases <- c("A","T","C","G")
    possibleBases <- foreach(i=1:length(klen)) %do% data.frame(gtools::permutations(n = length(bases), v = bases, r = klen[i], repeats.allowed = T))
    colsBases <- foreach(i=1:length(klen)) %do% colnames(possibleBases[[i]])
    foreach(i=1:length(klen)) %do% { possibleBases[[i]]$seq <- apply( possibleBases[[i]][ , colsBases[[i]] ] , 1 , paste , collapse = "" ) }
    foreach(i=1:length(klen)) %do% { residues[[i]] <- as.character(possibleBases[[i]]$seq) }
    
    #residues <- c("A","T","G","C")
    
    # Name the residues list by klen
    names(residues) <- klen
    
    # Remove unneeded vars
    rm(bases, possibleBases, colsBases)
    
    # List for each kPHMM Model constructed
    kphmmList <- list()
    
    # Message update
    cat(green(paste0("Starting model training...\n")))
    
    # kPHMM Model construction, iterate over each element of kMat
    foreach(i=1:length(kMat)) %do% {
      # residue vector is chosen depending on which klen is being used for the current model being constructed
      residue <- residues[[klenList[[i]]]]
      kphmmList[[i]] <- .kphmm_runTrainF(kMat[[i]], weights[[i]], residue, rankName[[i]],
                                         groupName[[i]], klenList[[i]], foldName[[i]], numMaxiter, 
                                         deltaLL, numCores)
    }
    
    # Double check a few plots to make sure the model output is good
    # plot(kphmmList[[1]]$model, main = paste0("Group: ", groupName[[1]], " Rank: ", rankName[[1]], sep = ""), arrexp = 1, textexp = 1)
    # plot(kphmmList[[2]]$model, main = paste0("Group: ", groupName[[2]], " Rank: ", rankName[[3]], sep = ""), arrexp = 1, textexp = 1)
    # plot(kphmmList[[3]]$model, main = paste0("Group: ", groupName[[3]], " Rank: ", rankName[[3]], sep = ""), arrexp = 1, textexp = 1)
    
    # Find which groups were successful in model construction and subset kMat & kphmmList by those groups
    modelTally <- unlist(sapply(kphmmList, function(x){ x$messages }))
    modelTally <- str_detect(modelTally, "^Success!")
    kphmmList <- kphmmList[modelTally]
    kMat <- kMat[modelTally]
    
    # remove unneeded vars
    rm(modelTally)
    
    # Create a df from vars extracted from kphmmList to indicate which groups were successful
    modelOverview <- list("Success" = data.frame("group" = sapply(kphmmList, function(x){ x$Group }),
                                                 "rank" = sapply(kphmmList, function(x){ x$Rank }),
                                                 "klen" = sapply(kphmmList, function(x){ x$klen }),
                                                 "fold" = sapply(kphmmList, function(x){ x$fold })),
                          "Fail" = kmerResults[['Tally']][which(!kmerResults[['Tally']]$group %in% sapply(kphmmList, function(x){ x$Group })),] %>% select(!numDiscrimKmers))
    
    # Writing models to files in the modelData folder
    
    # Creating directory for model files
    dirM <- paste0(dir, "/ModelData")
    dir.create(dirM, showWarnings = FALSE)
    
    # list for model names
    modelFiles <- list()
    # list for input matrices
    inputMatFiles <- list()
    
    # Iterate through each input matrix, writing each to disk
    foreach(i=1:length(kphmmList)) %do% {
      # Update progress  
      message('\r', green(paste0("Writing input matrix to .tsv file for group ", kphmmList[[i]]$Group, "                       ", sep="")), appendLF = FALSE)
      
      # Create filenames
      inputMatFiles[[i]] <- paste(dirM, "/input_matrix_", kphmmList[[i]]$Group, "_k", kphmmList[[i]]$klen, "_", kphmmList[[i]]$fold, ".tsv", sep="")
      
      # Also append filenames to each kphmmlist element
      kphmmList[[i]]$inputMat_filepath <- inputMatFiles[[i]]
      
      # write out each kphmm (encased in try statement if this function fails)
      try(write_tsv((as.data.frame(kMat[[i]]) %>% rownames_to_column()), col_names = TRUE, file = inputMatFiles[[i]]), silent = TRUE)
    }
    
    # Message update
    cat(green(paste0("All input matrices for each model have been written to the ", dirM, " directory in .tsv format.\n")))
    
    # Iterate through each model trained, writing each to disk
    foreach(i=1:length(kphmmList)) %do% {  
      # Update progress  
      message('\r', green(paste0("Writing model to .hmm file for group ", kphmmList[[i]]$Group, "                       ", sep="")), appendLF = FALSE)
      
      # Create filenames
      modelFiles[[i]] <- paste(dirM, "/model_", kphmmList[[i]]$Group, "_k", kphmmList[[i]]$klen, "_", kphmmList[[i]]$fold, ".hmm", sep="")
      
      # Also append filenames to each kphmmlist element
      kphmmList[[i]]$model_filepath <- modelFiles[[i]]
      
      # write out each kphmm (encased in try statement if this function fails)
      try(writePHMM(kphmmList[[i]]$model, file = modelFiles[[i]]), silent = TRUE)
    }
    
    # Remove unneeded vars
    rm(modelFiles, inputMatFiles)
    
    # Message update
    cat(green(paste0("All trained models have been written to the ", dirM, " directory in .hmm format.\n")))
    
    # Create a final list consisting of three elements: 
    # Model Overview - dataframe indicating which groups successfully had kPHMM models trained for them 
    # Model Data - Model kphmms for each group including the filepaths of each model and input matrix written to disk
    # Model Input Matrices - The input matrices used for each group during model setup
    modelList <- list("Overview" = modelOverview, "kphmmData" = kphmmList, "InputMatrices" = kMat)
    
    # Names for modelList (Model Data and Model Input Matrices)
    names(modelList[[2]]) <- paste0(sapply(kphmmList, function(x){ x$Group }), "_", 
                                    sapply(kphmmList, function(x){ x$klen }), "_", 
                                    sapply(kphmmList, function(x){ x$fold }), sep="")
    names(modelList[[3]]) <- paste0(sapply(kphmmList, function(x){ x$Group }), "_", 
                                    sapply(kphmmList, function(x){ x$klen }), "_", 
                                    sapply(kphmmList, function(x){ x$fold }), sep="")
    
    # Remove unneeded vars
    rm(kphmmList, kMat, modelOverview)
    
    # If cross-validation tests are NOT being performed
  } else {
    
    # Naming vectors for taxonomy and fold data
    nameVecTax <- c("recordID", ids)
    
    # List for kPosData
    kPosData <- list()
    
    # Grab the position data from the raw kmer data files
    foreach(i=1:length(klen)) %do% {
      # Read the data in chunkwise
      kPosData[[i]] <- read_chunkwise(as.character(kmerFiles[[i]]), sep="\t", chunk_size = chunkSize)
      # select position, kmers and taxonomy columns and collect from chunkwise objects
      kPosData[[i]] <- kPosData[[i]] %>%
        select(cols = c(3, 4, 6)) %>%
        # rename pos, kmers, taxonomy columns
        mutate(pos = cols1) %>%
        mutate(kmers = cols2) %>%
        mutate(taxonomy = cols3) %>%
        # select all relevant columns and collect data
        select(pos, kmers, taxonomy) %>%
        collect() %>% 
        separate(taxonomy, nameVecTax, sep = ";") %>%
        mutate(klen = klen[i])
    }
    
    # Rbind positional data across klengths
    kPosData <- data.frame(rbindlist(kPosData))
    
    # Divide the dataframe by subclassRank and reorganize columns for rank and group
    # Order
    kPosData_o <- kPosData %>% 
      select(pos, kmers, recordID, ordID, klen) %>% 
      mutate(group = ordID) %>% 
      mutate(rank = "Order") %>% 
      select(!ordID)
    
    # Family
    kPosData_f <- kPosData %>% 
      select(pos, kmers, recordID, famID, klen) %>% 
      mutate(group = famID) %>% 
      mutate(rank = "Family") %>% 
      select(!famID)
    
    # Genus
    kPosData_g <- kPosData %>% 
      select(pos, kmers, recordID, genID, klen) %>% 
      mutate(group = genID) %>% 
      mutate(rank = "Genus") %>% 
      select(!genID)
    
    # Species
    kPosData_s <- kPosData %>% 
      select(pos, kmers, recordID, spID, klen) %>% 
      mutate(group = spID) %>% 
      mutate(rank = "Species") %>% 
      select(!spID)
    
    # Combine back together with rbind
    kPosData <- rbind(kPosData_o, kPosData_f, kPosData_g, kPosData_s)
    
    # Remove unneeded vars
    rm(kPosData_o, kPosData_f, kPosData_g, kPosData_s)
    
    # Ensure both kmers and kPosData are in dataframe format
    kmers <- data.frame(kmers)
    kPosData <- data.frame(kPosData)
    
    # Compare groups here with groups in the kmers dataframe and subset out those that are not present (since those are unusable)
    kPosData <- kPosData[which(kPosData$group %in% kmers$group),]
    
    # Split into a list by folds, group & klen for both kPosData & kmers dataframes
    # kPosData dataframe
    kPosData <- kPosData %>% 
      group_by(group, klen) %>%
      group_split()
    
    # kmers dataframe
    kmers <- kmers %>%
      group_by(group, klen) %>%
      group_split()
    
    # This set of for loops will convert to DT for every group
    kmers <- foreach(i=1:length(kmers)) %do% data.table(kmers[[i]])
    kPosData <- foreach(i=1:length(kPosData)) %do% data.table(kPosData[[i]])
    
    # Setkeys for every dataframe across kPosData and kmers
    foreach(i=1:length(kmers)) %do% {
      setkey(kPosData[[i]], "kmers")
      setkey(kmers[[i]], "kmers")
    }
    
    # Merge kmer dataframes with kPosData dataframes using kmers as the key
    kPosData <- foreach(i=1:length(kPosData)) %do% merge.data.table(kPosData[[i]], kmers[[i]])
    
    # Find all unique kmers 
    uniqKmers <- foreach(i=1:length(kmers)) %do% unique(kmers[[i]]$kmers)
    
    # Count number of recordIDs
    recordNum <- foreach(i=1:length(kPosData)) %do% length(unique(kPosData[[i]]$recordID))
    
    # Tally up the number of disrim kmers at each position / total number of sequences
    kTally <- list()
    foreach(i=1:length(kPosData)) %do% {
      kTally[[i]] <- kPosData[[i]] %>%
        group_by(pos, kmers) %>%
        tally() %>% 
        mutate(n = n / recordNum[[i]]) 
    }
    
    # Search for positions with the highest numbers of discrim kmers, pick the top number of positions (as defined by numMaxKPos)
    foreach(i=1:length(kTally)) %do% {
      kTally[[i]] <- subset(kTally[[i]], kTally[[i]]$kmers %in% uniqKmers[[i]]) %>%
        arrange(desc(n)) %>%
        ungroup(pos, kmers) %>%
        slice_head(n = numMaxKPos)
    }
    
    # Subset by those positions (highest numbers of discrim kmers), 
    kMat <- list()
    foreach(i=1:length(kTally)) %do% {
      kMat[[i]] <- subset(kPosData[[i]], kPosData[[i]]$pos %in% kTally[[i]]$pos) %>%
        mutate(data = paste0(rank.x, ";", group.x, ";", klen.x, sep="")) %>%
        select(data, recordID, pos, kmers)
    }
    
    # Delete any duplicate rows (across all columns), relabel 'NC' record ids to be unique by rownumber and spread out data in a matrix like format for each list element, 
    # fill by gap symbol
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- dplyr::distinct(kMat[[i]]) %>%
        mutate(recordID = ifelse(recordID == "NC", paste0("NC", row_number(), sep=""), recordID)) %>%
        spread(pos, kmers, fill = "-")
    }
    
    # Find gappy kmer positions and filter them out, removal of these kmer positions is governed by the gapThreshold parameter
    kGap <- list()
    foreach(i=1:length(kMat)) %do% {
      kGap[[i]] <- kMat[[i]] %>%
        gather(pos, kmers) %>%
        filter(str_detect(pos,"^\\s*[0-9]*\\s*$")) %>%
        group_by(pos) %>% 
        dplyr::summarise(sumG = sum(kmers == "-") / nrow(kMat[[i]])) %>%
        filter(sumG <= gapThreshold)
    }
    
    # kmerCheck var to find which df's are empty from kGap
    kmerCheck <- foreach(i=1:length(kGap)) %do% {
      nrow(kGap[[i]]) == 0
    }
    
    # Subset kMat, kGap and uniqKmers by kmerCheck
    kMat <- kMat[!unlist(kmerCheck)]
    kGap <- kGap[!unlist(kmerCheck)]
    uniqKmers <- uniqKmers[!unlist(kmerCheck)]
    
    # Remove kmerCheck
    rm(kmerCheck)
    
    # Subset kMat by kGap positions
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- cbind.data.frame(kMat[[i]][,1:2], subset(kMat[[i]], select = colnames(kMat[[i]]) %in% kGap[[i]]$pos))
    }
    
    # Count how many discriminitive kmers each sequence has (identified by recordID) and filter by a minimum number of discrim kmers per recordID (as defined by numMinKSeq)
    kMat2 <- list()
    foreach(i=1:length(kMat)) %do% {
      kMat2[[i]] <- kMat[[i]] %>%
        mutate(count = across(.cols = 3:ncol(kMat[[i]]), ~ .x %in% uniqKmers[[i]])) %>%
        rowwise() %>%
        mutate(countDiscrim = across(.cols = contains("count"), .fns = sum)) %>%
        select(recordID, countDiscrim) %>%
        filter(countDiscrim > numMinKSeq[1])
    }
    
    # Add a weighted count column to assign weights to each sequence for model setup based on number of discrim kmer positions they posses
    # ex: a sequence with the all possible positions would be assigned the highest possible weight proportional to the total number of sequences
    # note: average weight across all sequences must equal 1 and sum of all weights must equal total number of sequences
    foreach(i=1:length(kMat2)) %do% {
      kMat2[[i]]$sumCount <- sum(kMat2[[i]]$countDiscrim)
      kMat2[[i]]$weight <- (kMat2[[i]]$countDiscrim / kMat2[[i]]$sumCount) * nrow(kMat2[[i]])
    }
    
    # Final subset by recordID to get the matrices that will be used in the PHHM's and add the previously calculated weight column
    foreach(i=1:length(kMat)) %do% {
      kMat[[i]] <- subset(kMat[[i]], kMat[[i]]$recordID %in% kMat2[[i]]$recordID)
      kMat[[i]]$weight <- kMat2[[i]]$weight
    }
    
    # Remove unneeded vars
    rm(kMat2, kGap, kTally, uniqKmers, recordNum, kPosData, kmers)
    
    # Transfer recordID column to rowname
    kMat <- foreach(i=1:length(kMat)) %do% {column_to_rownames(as.tibble(kMat[[i]]), var = "recordID")}
    
    # Convert kMat to matrix format
    foreach(i=1:length(kMat)) %do% {kMat[[i]] <- as.matrix(kMat[[i]])}
    
    # Check for groups that possess no discrim kmer positions (after the above steps) and filter them out
    kmerCheck <- foreach(i=1:length(kMat)) %do% {
      length(kMat[[i]]) == 0 | ncol(kMat[[i]]) < 3
    }
    
    # Subset kMat by kmerCheck
    kMat <- kMat[!unlist(kmerCheck)]
    
    # Remove unneeded vars
    rm(kmerCheck)
    
    # Take the first rowname element from each list element
    nameList <- unlist(foreach(i=1:length(kMat)) %do% {str_split(kMat[[i]][1,1], ";")}, recursive = FALSE)
    
    # rankName from nameList
    if(length(klen) != 1){
      rankName <- sapply(nameList, `[`, 1)
    } else {
      rankName <- as.list(sapply(nameList, `[`, 1))
    }
    
    # groupName from nameList
    if(length(klen) != 1){
      groupName <- sapply(nameList, `[`, 2)
    } else {
      groupName <- as.list(sapply(nameList, `[`, 2))
    }
    
    # klen from nameList
    if(length(klen) != 1){
      klenList <- sapply(nameList, `[`, 3)
    } else {
      klenList <- as.list(sapply(nameList, `[`, 3))
    }
    
    # Extract the weight column
    weights <- foreach(j=1:length(kMat)) %do% kMat[[j]][, ncol(kMat[[j]])]
    
    # Remove the data column and the weight column from each matrix 
    kMat <- foreach(i=1:length(kMat)) %do% kMat[[i]][,-1]
    kMat <- foreach(i=1:length(kMat)) %do% kMat[[i]][,-ncol(kMat[[i]])]
    
    # Create vectors of residues for the PHMM for each klen (all possible kmers for a given klen using gtools::permutations)
    residues <- list()
    bases <- c("A","T","C","G")
    possibleBases <- foreach(i=1:length(klen)) %do% data.frame(gtools::permutations(n = length(bases), v = bases, r = klen[i], repeats.allowed = T))
    colsBases <- foreach(i=1:length(klen)) %do% colnames(possibleBases[[i]])
    foreach(i=1:length(klen)) %do% { possibleBases[[i]]$seq <- apply( possibleBases[[i]][ , colsBases[[i]] ] , 1 , paste , collapse = "" ) }
    foreach(i=1:length(klen)) %do% { residues[[i]] <- as.character(possibleBases[[i]]$seq) }
    
    # Name the residues list by klen
    names(residues) <- klen
    
    # Remove unneeded vars
    rm(bases, possibleBases, colsBases)
    
    # List for each kPHMM Model constructed
    kphmmList <- list()
    
    # Message update
    cat(green(paste0("Starting model training...\n")))
    
    # kPHMM Model construction, iterate over each element of kMat
    foreach(i=1:length(kMat)) %do% {
      # residue vector is chosen depending on which klen is being used for the current model being constructed
      residue <- residues[[klenList[[i]]]]
      kphmmList[[i]] <- .kphmm_runTrainNF(kMat[[i]], weights[[i]], residue, rankName[[i]],
                                          groupName[[i]], klenList[[i]], numMaxiter, 
                                          deltaLL, numCores)
    }
    
    # Double check a few plots to make sure the model output is good
    # plot(kphmmList[[1]]$model, main = paste0("Group: ", groupName[[1]], " Rank: ", rankName[[1]], sep = ""), arrexp = 1, textexp = 1)
    # plot(kphmmList[[2]]$model, main = paste0("Group: ", groupName[[2]], " Rank: ", rankName[[3]], sep = ""), arrexp = 1, textexp = 1)
    # plot(kphmmList[[3]]$model, main = paste0("Group: ", groupName[[3]], " Rank: ", rankName[[3]], sep = ""), arrexp = 1, textexp = 1)
    
    # Find which groups were successful in model construction and subset kMat & kphmmList by those groups
    modelTally <- unlist(sapply(kphmmList, function(x){ x$messages }))
    modelTally <- str_detect(modelTally, "^Success!")
    kphmmList <- kphmmList[modelTally]
    kMat <- kMat[modelTally]
    
    # remove unneeded vars
    rm(modelTally)
    
    # Create a df from vars extracted from kphmmList to indicate which groups were successful
    modelOverview <- list("Success" = data.frame("group" = sapply(kphmmList, function(x){ x$Group }),
                                                 "rank" = sapply(kphmmList, function(x){ x$Rank }),
                                                 "klen" = sapply(kphmmList, function(x){ x$klen })),
                          "Fail" = kmerResults[['Tally']][which(!kmerResults[['Tally']]$group %in% sapply(kphmmList, function(x){ x$Group })),] %>% select(!numDiscrimKmers))
    
    # Writing models to files in the modelData folder
    
    # Creating directory for model files
    dirM <- paste0(dir, "/ModelData")
    dir.create(dirM, showWarnings = FALSE)
    
    # list for model names
    modelFiles <- list()
    # list for input matrices
    inputMatFiles <- list()
    
    # Iterate through each input matrix, writing each to disk
    foreach(i=1:length(kphmmList)) %do% {
      # Update progress  
      message('\r', green(paste0("Writing input matrix to .tsv file for group ", kphmmList[[i]]$Group, "                       ", sep="")), appendLF = FALSE)
      
      # Create filenames
      inputMatFiles[[i]] <- paste(dirM, "/input_matrix_", kphmmList[[i]]$Group, "_k", kphmmList[[i]]$klen, ".tsv", sep="")
      
      # Also append filenames to each kphmmlist element
      kphmmList[[i]]$inputMat_filepath <- inputMatFiles[[i]]
      
      # write out each kphmm (encased in try statement if this function fails)
      try(write_tsv((as.data.frame(kMat[[i]]) %>% rownames_to_column()), col_names = TRUE, file = inputMatFiles[[i]]), silent = TRUE)
    }
    
    # Message update
    cat(green(paste0("All input matrices for each model have been written to the ", dirM, " directory in .tsv format.\n")))
    
    # Iterate through each model that was trained, writing each to disk
    foreach(i=1:length(kphmmList)) %do% {  
      # Update progress  
      message('\r', green(paste0("Writing model to .hmm file for group ", kphmmList[[i]]$Group, "                       ", sep="")), appendLF = FALSE)
      
      # Create filenames
      modelFiles[[i]] <- paste(dirM, "/model_", kphmmList[[i]]$Group, "_k", kphmmList[[i]]$klen, ".hmm", sep="")
      
      # Also append filenames to each kphmmlist element
      kphmmList[[i]]$model_filepath <- modelFiles[[i]]
      
      # write out each kphmm (encased in try statement if this function fails)
      try(writePHMM(kphmmList[[i]]$model, file = modelFiles[[i]]), silent = TRUE)
    }
    
    # Remove unneeded vars
    rm(modelFiles, inputMatFiles)
    
    # Message update
    cat(green(paste0("All trained models have been written to the ", dirM, " directory in .hmm format.\n")))
    
    # Create a final list consisting of three elements: 
    # Model Overview - dataframe indicating which groups successfully had kPHMM models trained for them 
    # Model Data - Model kphmms for each group including the filepaths of each model and input matrix written to disk
    # Model Input Matrices - The input matrices used for each group during model setup
    modelList <- list("Overview" = modelOverview, "kphmmData" = kphmmList, "InputMatrices" = kMat)
    
    # Names for modelList (Model Data and Model Input Matrices)
    names(modelList[[2]]) <- paste0(sapply(kphmmList, function(x){ x$Group }), "_", 
                                    sapply(kphmmList, function(x){ x$klen }), sep="")
    names(modelList[[3]]) <- paste0(sapply(kphmmList, function(x){ x$Group }), "_", 
                                    sapply(kphmmList, function(x){ x$klen }), sep="")
    
    # Remove unneeded vars
    rm(kphmmList, kMat, modelOverview)
  }
  return(modelList)
}

kphmm_ClassifyCV <- function(kmerFiles, modelList, numFolds, subclassRank, klen, classRank){
  # Stop the function if numFolds is less than 3 (since that would mean cv tests were not being performed)
  if (numFolds < 3) {
    stop("Ensure that the numFolds numeric is greater or equal to 3 to use this function. To run classifications without performing cross-validation tests,
          please use kphmm_ClassifyI (PHMM_Test.R file) instead.")
  }
  
  # Vector for subclassRank column ids
  ids <- paste0(tolower(ifelse(subclassRank != "Species", substr(subclassRank, 1, 3), substr(subclassRank, 1, 2))), "ID", sep="")
  
  # Naming vectors for taxonomy and fold data
  nameVecFolds <- paste0("f", 1:numFolds, sep="")
  nameVecTax <- c("recordID", ids, nameVecFolds)
  
  # List for kPosData
  kPosData <- list()
  
  # Grab the position data from the raw kmer data files
  foreach(i=1:length(klen)) %do% {
    # Read the data in chunkwise
    kPosData[[i]] <- read_chunkwise(as.character(kmerFiles[[i]]), sep="\t", chunk_size = chunkSize)
    
    # select position, kmers and taxonomy columns and collect from chunkwise objects
    kPosData[[i]] <- kPosData[[i]] %>%
      select(cols = c(3, 4, 6)) %>%
      # rename pos, kmers, taxonomy columns
      mutate(pos = cols1) %>%
      mutate(kmers = cols2) %>%
      mutate(taxonomy = cols3) %>%
      # select all relevant columns and collect data
      select(pos, kmers, taxonomy) %>%
      collect() %>% 
      separate(taxonomy, nameVecTax, sep = ";") %>%
      mutate(klen = klen[i])
  }
  
  # Rbind positional data across klengths
  kPosData <- data.frame(rbindlist(kPosData))
  
  # list for folds
  folds <- replicate(numFolds, list())
  
  # Find the rows labelled as "test" from each fold
  foreach(i=1:numFolds) %do% {
    # Subset the data into different folds
    folds[[i]] <- kPosData[which(kPosData[,i+7] == "test"),]
    # modify kmers column to also include fold number
    folds[[i]]$folds <- paste0("f", i)
  }
  
  # Rbind all folds together into one large dataframe and remove duplicates
  kPosData <- rbindlist(folds) %>% 
    select(folds, pos, kmers, recordID, ordID, famID, genID, spID, klen) %>%
    distinct()
  
  # Divide the dataframe by subclassRank and reorganize columns for rank and group
  # Order
  kPosData_o <- kPosData %>% 
    select(folds, pos, kmers, recordID, ordID, klen) %>% 
    mutate(group = ordID) %>% 
    mutate(rank = "Order") %>% 
    select(!ordID)
  
  # Family
  kPosData_f <- kPosData %>% 
    select(folds, pos, kmers, recordID, famID, klen) %>% 
    mutate(group = famID) %>% 
    mutate(rank = "Family") %>% 
    select(!famID)
  
  # Genus
  kPosData_g <- kPosData %>% 
    select(folds, pos, kmers, recordID, genID, klen) %>% 
    mutate(group = genID) %>% 
    mutate(rank = "Genus") %>% 
    select(!genID)
  
  # Species
  kPosData_s <- kPosData %>% 
    select(folds, pos, kmers, recordID, spID, klen) %>% 
    mutate(group = spID) %>% 
    mutate(rank = "Species") %>% 
    select(!spID)
  
  # Combine back together with rbind
  kPosData <- rbind(kPosData_o, kPosData_f, kPosData_g, kPosData_s)
  
  # Remove unneeded vars
  rm(folds, kPosData_o, kPosData_f, kPosData_g, kPosData_s)
  
  # Create a dataframe of all unique record ids along with their taxonomy that can be cross-referenced to verify taxonomic classifications
  recordDict <- kPosData %>% 
    filter(recordID != "NC") %>% 
    select(group, rank, recordID) %>% 
    distinct()
  
  # Dataframe divided by fold which will be used for the testing of all sequences labelled "test" on trained kphmm models
  kTest <- kPosData %>% 
    filter(recordID != "NC") %>%
    group_by(folds) %>%
    group_split()
  
  # Take dataframe divided by fold and subdivide again by klen
  kTest <- foreach(i=1:length(kTest)) %do% {
    kTest[[i]] %>% 
      group_by(klen) %>%
      group_split()
  }
  
  # Subdivide again by rank
  kTest <- foreach(i=1:length(kTest)) %do% {
    foreach(j=1:length(kTest[[i]])) %do% {
      kTest[[i]][[j]] %>% 
        group_by(rank) %>%
        group_split()
    }
  }
  
  # Finally subdivide by recordID, eliminate duplicate rows and arrange kmers by position, select by recordID, kmers and pos and then name different levels of the list
  foreach(i=1:length(kTest)) %do% {
    foreach(j=1:length(kTest[[i]])) %do% {
      kTest[[i]][[j]] <- foreach(l=1:length(kTest[[i]][[j]])) %do% {
        kTest[[i]][[j]][[l]] %>% 
          group_by(recordID) %>%
          distinct() %>% 
          arrange(pos) %>%
          select(recordID, kmers, pos) %>%
          group_split()
      }
      names(kTest[[i]][[j]]) <- sort(subclassRank)
    }
    names(kTest[[i]]) <- paste0(klen[1:length(klen)], sep="")
  }
  names(kTest) <- paste0("f", 1:numFolds, sep="")
  
  # Add rank to names of kphmm's and kMat's
  ranks <- unlist(foreach(i=1:length(modelList$kphmmData)) %do% modelList$kphmmData[[i]]$Rank)
  names(modelList$kphmmData) <- paste0(names(modelList$kphmmData), "_", ranks, sep="")
  names(modelList$InputMatrices) <- paste0(names(modelList$InputMatrices), "_", ranks, sep="")
  
  # remove unneeded vars
  rm(kPosData)
  
  # Organize kphmm data the same way as with the test data
  modelTest <- replicate(numFolds, list())
  subclassRankT <- sort(subclassRank)
  foreach(i=1:length(modelTest)) %do% {
    modelTest[[i]] <- replicate(length(klen), list())
    foreach(j=1:length(klen)) %do% {
      modelTest[[i]][[j]] <- replicate(length(subclassRank), list())
      foreach(l=1:length(subclassRank)) %do% {
        modelTest[[i]][[j]][[l]] <- modelList$kphmmData[str_detect(names(modelList$kphmmData), paste0("_", klen[j], "_f", i, "_", subclassRankT[l], sep=""))]
      }
      names(modelTest[[i]][[j]]) <- subclassRankT
    }
    names(modelTest[[i]]) <- paste0(klen[1:length(klen)], sep="")
  }
  names(modelTest) <- paste0("f", 1:numFolds, sep="")
  
  # Organize matrix data the same way it is organized for the test data and model test data
  matTest <- replicate(numFolds, list())
  foreach(i=1:length(matTest)) %do% {
    matTest[[i]] <- replicate(length(klen), list())
    foreach(j=1:length(klen)) %do% {
      matTest[[i]][[j]] <- replicate(length(subclassRank), list())
      foreach(l=1:length(subclassRank)) %do% {
        matTest[[i]][[j]][[l]] <- modelList$InputMatrices[str_detect(names(modelList$InputMatrices), paste0("_", klen[j], "_f", i, "_", subclassRankT[l], sep=""))]
      }
      names(matTest[[i]][[j]]) <- subclassRankT
    }
    names(matTest[[i]]) <- paste0(klen[1:length(klen)], sep="")
  }
  names(matTest) <- paste0("f", 1:numFolds, sep="")
  
  # Test each test sequence using the models from modelTest and the matrices from matTest
  # Name each sublist with its corresponding name (either fold, klen, subclassrank, recordID or modelName)
  modelResult <- replicate(numFolds, list())
  foreach(i=1:length(modelResult)) %do% {
    modelResult[[i]] <- replicate(length(klen), list())
    # Subdivide by klen
    foreach(j=1:length(klen)) %do% {
      modelResult[[i]][[j]] <- replicate(length(subclassRank), list())
      # Subdivide by subclassRank
      foreach(l=1:length(subclassRank)) %do% {
        #*** Added an if statement here in case matTest is empty
        if(length(matTest[[i]][[j]][[l]]) > 0){
          modelResult[[i]][[j]][[l]] <- replicate(length(kTest[[i]][[j]][[l]]), list())
          # Subdivide by recordID
          foreach(x=1:length(modelResult[[i]][[j]][[l]])) %do% {
            # Update progress  
            message('\r', green(paste0("Performing classifications for test sequence: ", kTest[[i]][[j]][[l]][[x]]$recordID[1], "                       ", sep="")), appendLF = FALSE)
            # use the classify_ForwardCV function on each recordID testing it with each matrix and model and return the top scoring result in a dataframe
            modelResult[[i]][[j]][[l]][[x]] <- .kphmm_classifyForwardCV(kTest[[i]][[j]][[l]][[x]], 
                                                                        matTest[[i]][[j]][[l]], 
                                                                        modelTest[[i]][[j]][[l]],
                                                                        kTest[[i]][[j]][[l]][[x]]$recordID[1])
          }
        } else {
          modelResult[[i]][[j]][[l]] <- NULL
        }
      }
    }
  }
  
  # Rbind each list of each modelResult dataframe, resulting in one single dataframe
  modelResult <- foreach(i=1:numFolds) %do% {
    modelResult[[i]] <- foreach(j=1:length(klen)) %do% {
      modelResult[[i]][[j]] <- foreach(l=1:length(subclassRank)) %do% {
        rbindlist(modelResult[[i]][[j]][[l]])
      }
      rbindlist(modelResult[[i]][[j]])
    }
    rbindlist(modelResult[[i]])
  }
  modelResult <- rbindlist(modelResult)
  
  # setkeys (recordID) before merging
  setkey(modelResult, "rank", "recordID")
  setkey(recordDict, "rank", "recordID")
  # Merge together
  modelResult <- merge.data.table(modelResult, recordDict)
  # Column names
  colnames(modelResult)[4] <- "assignedTaxa"
  colnames(modelResult)[7] <- "Group"
  
  # Use the model results to determine percentage correct taxonomic assignments on a per group level
  percentCorrect_group <- suppressWarnings(modelResult %>%
                                           group_by(rank, Group, fold, klen) %>%
                                           dplyr::summarize(`CorrectId%` = (sum(assignedTaxa == Group) / n()) * 100, n = n()) %>% 
                                           ungroup(rank, Group, fold, klen) %>%
                                           mutate(`Fold&klength` = paste0(fold, ",", klen, sep="")))
  # install.packages("ROCR")
  # library(ROCR)
  
  # Convert to probability score and then using a user-specified conf cutoff find TP, FP, TN, FN values
  confidenceVal_group <- list()
  confVal <- seq(0.01, 1, by=0.01)
  foreach(i=1:100) %do% {
    confidenceVal_group[[i]] <- modelResult %>%
      mutate(probScore = exp(score)/(1+exp(score))) %>%
      mutate(TP = probScore >= confVal[i] & assignedTaxa == Group) %>%
      mutate(FP = probScore >= confVal[i] & assignedTaxa != Group) %>%
      mutate(TN = probScore < confVal[i] & assignedTaxa != Group) %>%
      mutate(FN = probScore < confVal[i] & assignedTaxa == Group) %>%
      mutate("confThreshold" = confVal[i])
  }
  
  # Tally up TPs across conf thresholds
  tpTally <- list()
  foreach(i=1:100) %do% {
    tpTally[[i]] <- confidenceVal_group[[i]] %>%
      group_by(rank, klen, fold) %>%
      dplyr::summarize("nTruePositive" = sum(TP == TRUE), "nTotal" = n()) %>%
      mutate("%TP" = (nTruePositive / nTotal) * 100) %>%
      mutate("confThreshold" = confVal[i])
  }
  
  # Tally up FPs across conf thresholds
  fpTally <- list()
  foreach(i=1:100) %do% {
    fpTally[[i]] <- confidenceVal_group[[i]] %>%
      group_by(rank, klen, fold) %>%
      dplyr::summarize("nFalsePositive" = sum(FP == TRUE), "nTotal" = n()) %>%
      mutate("%FP" = (nFalsePositive / nTotal) * 100) %>%
      mutate("confThreshold" = confVal[i])
  }
  
  # Tally up FPs across conf thresholds
  tnTally <- list()
  foreach(i=1:100) %do% {
    tnTally[[i]] <- confidenceVal_group[[i]] %>%
      group_by(rank, klen, fold) %>%
      dplyr::summarize("nTrueNegative" = sum(TN == TRUE), "nTotal" = n()) %>%
      mutate("%TN" = (nTrueNegative / nTotal) * 100) %>%
      mutate("confThreshold" = confVal[i])
  }
  
  # Tally up FPs across conf thresholds
  fnTally <- list()
  foreach(i=1:100) %do% {
    fnTally[[i]] <- confidenceVal_group[[i]] %>%
      group_by(rank, klen, fold) %>%
      dplyr::summarize("nFalseNegative" = sum(FN == TRUE), "nTotal" = n()) %>%
      mutate("%FN" = (nFalseNegative / nTotal) * 100) %>%
      mutate("confThreshold" = confVal[i])
  }
  
  # Tally up % correct ids
  percentTally <- modelResult %>%
    group_by(rank, klen, fold) %>%
    dplyr::summarize("correctID" = sum(assignedTaxa == Group), "nTotal" = n()) %>%
    mutate("%CorrectID" = (correctID / nTotal) * 100)
  
  # tally's combined for FP, FN, TN, TP across confidence thresholds
  tallyCombined <- list()
  foreach(i=1:100) %do% {
    tallyCombined[[i]] <- cbind(tpTally[[i]][,c(1,2,3,5,7)], tpTally[[i]][,c(4,6)], fpTally[[i]][,c(4,6)], tnTally[[i]][,c(4,6)], fnTally[[i]][,c(4,6)])
  }
  
  # rbind tallyCombined and confidenceVal_group
  tallyCombined <- rbindlist(tallyCombined)
  confidenceVal_group <- rbindlist(confidenceVal_group)
  
  # Create a final list of the results with all of the above dataframes
  modelResults <- list("ResultsData" = list("CorrectID% - Rank" = percentTally,
                                            "TP,FP,TN,FN% - Rank" = tallyCombined,
                                            "CorrectID% - Group" = percentCorrect_group,
                                            "TP,FP,TN,FN - Group" = confidenceVal_group),
                       "PlotFilepaths" = list("Order - HeatMap" = NULL,
                                              "Family - HeatMap" = NULL,
                                              "Genus - HeatMap" = NULL,
                                              "Species - HeatMap" = NULL,
                                              "Order - ROC,MCC,F1" = NULL,
                                              "Family - ROC,MCC,F1" = NULL,
                                              "Genus - ROC,MCC,F1" = NULL,
                                              "Species - ROC,MCC,F1" = NULL))
  
  # Return the model results
  return(modelResults)
}

# Internal function to subset by position according to each input matrix and use the forward algorithm on each kphmm with each set of kmer data (corresponding to one recordID each)
.kphmm_classifyForwardCV <- function(kTest, matTest, modelTest, recordID){

  # Subset positions by matrix positions
  kTestSub <- foreach(i=1:length(matTest)) %do% { subset(kTest, kTest$pos %in% colnames(matTest[[i]])) }
  
  # Extract kmers and name them with positions
  kTestSubKmers <- foreach(i=1:length(modelTest)) %do% { as.character(kTestSub[[i]]$kmers) }
  foreach(i=1:length(modelTest)) %do% { names(kTestSubKmers[[i]]) <- as.numeric(kTestSub[[i]]$pos) }
  
  # Run the forward algorithm to test the subsetted test kmer data with each model
  kTestProb <- foreach(i=1:length(modelTest)) %do% { forward(modelTest[[i]]$model, kTestSubKmers[[i]], odds = TRUE) }
  
  # Convert results data to dataframe format
  kTestProb <- foreach(i=1:length(modelTest)) %do% { data.frame("recordID" = recordID,
                                                                "score" = kTestProb[[i]]$score, 
                                                                "group" = modelTest[[i]]$Group, 
                                                                "rank" = modelTest[[i]]$Rank, 
                                                                "klen" = modelTest[[i]]$klen, 
                                                                "fold" = modelTest[[i]]$fold) }
  # rbind the list of dataframes
  kTestProb <- rbindlist(kTestProb)
  
  # Select the top scoring groups (the number of which is determined by the numTopGroups param, default = 1)  
  kTestProb <- kTestProb %>% slice_max(score, n = 1)

  # Return kTestProb
  return(kTestProb)
}

# Main function to classify sequences independently from sequences provided by the user on trained kphmm models (which may or may have not been trained by the user)
# modelResults <- kphmm_ClassifyI(subclassRank){
# }
