#!/usr/bin/env Rscript
library(minfi)

start.time <- Sys.time()

# get command line input
loc.out <- "/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Data"
loc.data <- "/groups/umcg-griac/tmp03/rawdata/450k/Maki_850K/ImageData"

# loading data
targets <- read.metharray.sheet(file.path("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Samplesheets"), pattern = "*.csv", recursive = TRUE)
# fixing an issue where the basenames are not specified
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})
# adding the sub directories too the paths with file.path basename = targets Slide, recursive is not enough
rg.set <- read.metharray.exp(base=file.path(loc.data, basename=targets$Slide), targets=targets, recursive=TRUE)
#saving the rg set to load it later
save(rg.set, targets, file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
