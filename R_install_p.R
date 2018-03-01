#!/usr/bin/env Rscript
# installation of bioconductor, minfi and loading of libraries
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("minfi")
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

