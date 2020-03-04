#!/usr/bin/env Rscript
install.packages(pkgs="plyr",repos = "http://cran.us.r-project.org")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
install.packages("tools/VirFinder/VirFinder_1.1.tar.gz", repos = NULL, type="source")
library(VirFinder)

fasta <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

predResult <- VF.pred(fasta)
# predResult$qvalue <- VF.qvalue(predResult$pvalue)
write.table(predResult, out, sep='\t', row.names=F, quote=F)