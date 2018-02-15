#!/usr/bin/env Rscript

library(minfi)

# load in cml arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are at least two arguments: if not, return an error
if (length(args)!= 2) {
  stop("This script needs an output dir and an input dir.n", call.=FALSE)
}

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log_0.txt")
start.time <- Sys.time()

# get command line input
loc.out <- args[1]
loc.data <- args[2]

# loading data
targets <- read.metharray.sheet(file.path(paste(loc.data, "SampleSheets", sep = "/")))
rg.set <- read.metharray.exp(file.path(paste(loc.data, "ImageData", sep = "/"), targets=targets))
save(rg.set, targets, file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
