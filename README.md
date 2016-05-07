# Co-expression Network Tool

A tool that calculates co-expressed gene network from GCT file of gene expressions.

For more complete documentation, see [the Yale CBB752 course website](http://cbb752spring2016.github.io/Network). 

Note: This tool is part of a set of bioinformatic and biological structure tools created for CBB752 at Yale University in the Spring 2016. The website containing links to the set of tools can be found at: <https://github.com/CBB752Spring2016/CBB752Spring2016.github.io>.

***


### Usage:


```
Rscript GCT_to_net.R -i (input .txt file in GCT format) -p (p value in [0,1]) -o (output .csv file name)
```

Example:    

```
Rscript GCT_to_net.R -i input.gct -p 0.01 -o output.csv
```

Note:       

* Input: GCT file, which is a text file with 1st row = column titles
              
* Columns: Name, Description, Experiment1, Experiment 2, exp3, exp4 ...
              
* One row per gene
              
* Output: a csv. of pairwise gene interactions, Pearson R values, p-values and Bonferroni-adjusted p-values, comma-delimited

* example.gct: constructed from all_aml_train.gct, remove 2 header lines, keep first 100 genes

***

This tool has also been implemented [in python by EdKong](https://github.com/EdKong/CBB752_Final_Project_3.1).

