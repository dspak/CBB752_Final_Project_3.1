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

# Load the required packages
#install.packages("argparser")
library(argparser)
if(!require(Hmisc)){install.packages("Hmisc")}


# Set to run as Rscript GCT_to_net.R /path/to/inputfile.gct
# To run interactively, use the setwd and input files below
args = commandArgs(trailingOnly=TRUE)
argv <- parse_args(p)


# specify our desired options
# by default ArgumentParser will add an help option
i <- add_argument(i, "--input", default=TRUE,
                    help='input file, tab delimited in a .txt file, MI TAB 2.5 format')
parser$add_argument("-p", "--pvalue", action="store_false", default=TRUE,
                    dest="verbose", help='input pvalue between 0 and 1')
parser$add_argument("-o", "--output", action="store_false",
                    dest="verbose", help='output filename')

input <- read.csv(file = args[1], header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F)


# To run interactively, uncomment the following two lines
setwd("~/Box Sync/coursework/CBB752_BioinformaticsMiningSimulation/final/CBB752_Final_Project_3.1/")
input <- read.table(file = "all_aml_train.preprocessed.gct", header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F, row.names = 1)




### Excluded from final product, include if using raw (not preprocessed) data
# ~~~~~~~~~~~~~
# ~~ STEP 1) Apply low and high thresholds and remove values that never change
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


# ~~~~~~~~~~~~~
# ~~ STEP 2) Find co-expressed genes by Pearson correlation
# ~~~~~~~~~~~~~



n <- ncol(input)-2

# http://stats.stackexchange.com/questions/153937/finding-p-value-in-pearson-correlation-in-r

r <- rcorr(t(input[,3:ncol(input)]), type = "pearson")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
df <- flattenCorrMatrix(r$r, r$P)
df$bonferroni <- p.adjust(df$p, method = "bonferroni")
df$t.statistic <- df$cor / (sqrt( (1-(df$cor^2)/(ncol(input)-4) )))

output <- subset(df, bonferroni <= 0.05)
output2 <- data.frame(output$row, output$column, output$cor, output$t.statistic, output$p, output$bonferroni)
colnames(output2) <- c("InteractorA", "InteractorB", "Correlation","t-statistic", "p-value", "Bonferroni-adjusted p-value")
write.csv(output2, "output_coexpressed_R.csv", row.names = F, quote = F, col.names = T)










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
# 
# 
# 
