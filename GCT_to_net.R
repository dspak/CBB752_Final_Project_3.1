#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title :  GCT_to_net.R 
# Version : 1.0
#
# Purpose : A tool that calculates co-expressed gene network from GCT file
#           of gene expressions. Uses a Pearson correlation and then a t-test
#           to prioritize the correlations
#  
# Version Notes : 
#
# Created.date  : 27 Apr 2016
# Created.by    : Dan Spakowicz
# Updated.date  :  
# Updated.by    : 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Usage:      Rscript GCT_to_net.R /path/to/inputfile.gct /path/to/outputfile.csv
# Example:    python GCT_to_net.py -i input.gct -p 0.01 -o output.csv
# Note:       Input: GCT file, which is a text file with 1st row = column titles
#             Columns: Name, Description, Experiment1, Experiment 2, exp3, exp4 ...
#             One row per gene
#             Output: a csv. of pairwise gene interactions, Pearson R values, and p-values, comma-delimited
#
#             example.gct: constructed from all_aml_train.gct, remove 2 header lines, keep first 10 genes


# Set to run as Rscript GCT_to_net.R /path/to/inputfile.gct
# To run interactively, use the setwd and input files below
args = commandArgs(trailingOnly=TRUE)

input <- read.csv(file = args[1], header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F)


# To run interactively, uncomment the following two lines
setwd("~/Box Sync/coursework/CBB752_BioinformaticsMiningSimulation/final/CBB752_Final_Project_3.1/")
input <- read.table(file = "all_aml_train.gct", header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F)

# Load the required packages


# ~~~~~~~~~~~~~
# ~~ STEP 1) Apply low and high thresholds and remove values that never change
# ~~~~~~~~~~~~~
# This requires first cleaning the data. As done using GenePattern earlier in this assignment, values <20 will be set to 20, values > 20000 will be set to 20000, and values that change by less than 3 fold between any sample will be removed.

# set values lower than 20 to 20
LowThreshold <- function(x) { x[x<20] <- 20; x }
df.low <- as.data.frame(LowThreshold(input))

#set values higher than 20000 to 20000
HighThreshold <- function(x) { x[x>20000] <- 20000; x }
df.low.high <- as.data.frame(HighThreshold(df.low[3:40]))

# Remove rows that vary by less than 3 fold in any sample
df.low.high$folddiff <- apply(df.low.high, 1, FUN = function(x) {max(x)/min(x)})
df.desc <- cbind(df.low$Name, df.low$Description, df.low.high)
df.clean <- df.desc[df.low.high$folddiff > 3,]
df.clean$folddiff <- NULL


# ~~~~~~~~~~~~~
# ~~ STEP 2) Find co-expressed genes by Pearson correlation
# ~~~~~~~~~~~~~

p.cor <- cor(df.clean[,3:ncol(df.clean)], method = "pearson")



