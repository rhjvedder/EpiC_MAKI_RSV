#!/usr/bin/env Rscript
source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
library(minfi)

start.time <- Sys.time()

# get command line input
loc.out <- "/data/f114798/Data"
loc.data <- "/data/f114798/Data/ImageData"

# loading data
targets <- read.metharray.sheet(file.path("/data/f114798/Samplesheets"), pattern = "*.csv", recursive = TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})
rg.set <- read.metharray.exp(base=file.path(loc.data, basename=targets$Slide), targets=targets, recursive=TRUE)
save(rg.set, targets, file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
