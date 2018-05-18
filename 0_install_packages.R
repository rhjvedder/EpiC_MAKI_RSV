#!/usr/bin/env Rscript
# installation of bioconductor, minfi and loading of libraries
source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
biocLite("GenABEL")
biocLite("WGCNA")
install.packages(c("ggplots2", "ggfortify", "RColorBrewer", "gridExtra", "pheatmap", "sva", "MASS", "compare", "tableone", "matrixStats", "plyr", "qqman", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))