## A strobemer-Based Profile Hidden Markov Model (sPHMM) approach to Taxonomic Identification of DNA Barcode Data.

#### To run this pipeline, the following files are needed in the same working directory in R/RStudio (R version 4.4.1 or higher):
sPHMM_Pipeline_Train&CV.R, sPHMM_Functions.R, strobemer_extract.cpp, strobemer_chop.cpp, 
strobemer_extract.h, strobemer_filter.cpp, strobemer_filter_unique.cpp and must also include a 
folder containing the fasta files to be run on this pipeline ex: seqData/

#### Functions get loaded from a separate R script and must be loaded before the analysis proceeds
source("sPHMM_Functions.R")


