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
# Updated.date  : 29 Apr 2016 
# Updated.by    : DS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Usage:      Rscript GCT_to_net.R --input /path/to/inputfile.gct --pvalue 0.05 -output /path/to/outputfile.csv
# Example:    python GCT_to_net.py -i input.gct -p 0.05 -o output.csv
# Note:       Input: GCT file, which is a text file with 1st row = column titles
#             Columns: Name, Description, Experiment1, Experiment 2, exp3, exp4 ...
#             One row per gene
#             Output: a csv. of pairwise gene interactions, Pearson R values, and p-values, comma-delimited
#
#             example.gct: constructed from all_aml_train.gct, remove 2 header lines, keep first 10 genes

# Load the required packages

if(!require(optparse)){install.packages("optparse", 
                                        repos = "http://cran.us.r-project.org")}
require(optparse)
if(!require(Hmisc)){install.packages("Hmisc", 
                                     repos = "http://cran.us.r-project.org")}
require(Hmisc)

# set arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--pvalue"), type="character", default=0.05, 
              help="input pvalue between 0 and 1", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              default="output_coexpressed_R.csv", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# To run interactively, uncomment the following two lines
# setwd("~/Box Sync/coursework/CBB752_BioinformaticsMiningSimulation/final/CBB752_Final_Project_3.1/")
# infile <- read.table(file = "all_aml_train.preprocessed.gct", header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F, row.names = 1)

# ~~~~~~~~~~~~~
# ~~ Find co-expressed genes by Pearson correlation, output table
# ~~~~~~~~~~~~~

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame (
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = cormat[ut],
    p = pmat[ut])
}

PearsonCorrelationTable <- function(infile, pvalue, outfile) {
  
  # read in data matrix
  infile <- read.csv(file = opt$input, header = T, sep = "\t", 
                     strip.white = T, quote = "\"", stringsAsFactors = F)
  
  # pearson correlation excluding name columns using Hmisc
  r <- rcorr(t(infile[,3:ncol(infile)]), type = "pearson")
  
  # reformat matrix to two columns of interactions with values
  df <- flattenCorrMatrix(r$r, r$P)
  
  # correct for multiple hypothesis testing
  df$bonferroni <- p.adjust(df$p, method = "bonferroni")
  
  # calculate the t-statistic to include in the final output
  df$t.statistic <- df$cor / (sqrt( (1-(df$cor^2)/(ncol(infile)-4) )))
  
  # subset the result to only significant correlations
  significants <- subset(df, bonferroni <= pvalue)
  
  # reformat output
  output <- data.frame(significants$row, significants$column, significants$cor, 
                        significants$t.statistic, significants$p, 
                        significants$bonferroni)
  colnames(output) <- c("InteractorA", "InteractorB", "Correlation",
                        "t-statistic", "p-value", "Bonferroni-adjusted p-value")
  
  # write to file
  write.csv(output, outfile, 
            row.names = F, quote = F, col.names = T)
}

# function to reformat matrix to two columns of interactions with values



PearsonCorrelationTable(opt$input, opt$pvalue, opt$out)





# http://stats.stackexchange.com/questions/153937/finding-p-value-in-pearson-correlation-in-r

#Excluded from final version, include if using unpreprocessed data
# ~~~~~~~~~~~~~ 
# ~~ Apply low and high thresholds and remove values that never change
# ~~~~~~~~~~~~~
# This requires first cleaning the data. As done using GenePattern earlier in this assignment, values <20 will be set to 20, values > 20000 will be set to 20000, and values that change by less than 3 fold between any sample will be removed.
# 
# # set values lower than 20 to 20
# LowThreshold <- function(x) { x[x<20] <- 20; x }
# df.low <- as.data.frame(LowThreshold(input))
# 
# #set values higher than 20000 to 20000
# HighThreshold <- function(x) { x[x>20000] <- 20000; x }
# df.low.high <- as.data.frame(HighThreshold(df.low[3:40]))
# 
# # Remove rows that vary by less than 3 fold in any sample
# df.low.high$folddiff <- apply(df.low.high, 1, FUN = function(x) {max(x)/min(x)})
# df.desc <- cbind(df.low$Name, df.low$Description, df.low.high)
# df.clean <- df.desc[df.low.high$folddiff > 3,]
# df.clean$folddiff <- NULL




# 
# names(r$r)
#   tstat = r * sqrt((n-2)/(1-r*r))
#   pval = pt(abs(tstat), df = n-2, lower = FALSE)*2
#   head(input)
#   colnames(pval) <- input$Description
# mtcars
#   pval2 <- data.frame(input$Description, pval)
#   melt <- melt(pval2)
#   head(melt)
#   
#   melt$adjusted <- p.adjust(melt$value, method = "bonferroni")
#   hist(melt$adjusted)
#   
#   table(melt$adjusted < 0.05)['TRUE'] / table(melt$adjusted < 0.05)['FALSE']
#   
# head (pval)
# head(FWER)
# hist(out)
# 
# # calculate t-stat, n-2 degrees of freedom
#   tstat = r*numpy.sqrt((n-2)/(1-r*r))
# 
# # find p-value for the double-sided test. Students t, n-2 degrees of freedom
# pval = stats.t.sf(numpy.abs(tstat), n-2)*2
