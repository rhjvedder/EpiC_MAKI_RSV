#!/usr/bin/env Rscript
library(minfi)

# This file contains functions for the pipeline to process 850k epic illumina data
# loading the rg.set
loc.idat <- "/data/f114798/Data/ImageData"
loc.sheet <- "/data/f114798/Samplesheets"
loc.out <- "/data/f114798/Data"

start.time <- Sys.time()

# conversion and storing of idat to rg files with option to subset on bad samples
targets <- read.metharray.sheet(loc.sheet, pattern="*.csv", recursive=TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
det.p <- detectionP(rg.set)
save(det.p, file=paste(loc.out, "Detection_Pvalues.Rdata", sep="/"))
# remove sample if more then 20% of probes failed on the sample
bad.samples <- colMeans(det.p > 0.01) > 0.05
bad.sample.names.detP <- colnames(det.p[,bad.samples])
bad.sample.detp.rs <- data.frame(CallRate=colMeans(det.p[,bad.samples] > 0.01))
rownames(bad.sample.detp.rs) <- rg.set$Sample_Name[rg.set$Basename %in% bad.sample.names.detP]
write.csv(bad.sample.detp.rs, file=paste(loc.out, "Bad_Samples_CallRate.csv", sep="/"))
targets <- targets[!bad.samples,]
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(rg.set, targets, bad.sample.names.detP, bad.samples, file=paste(loc.out, "RG_Channel_Set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
