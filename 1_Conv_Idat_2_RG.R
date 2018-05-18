#!/usr/bin/env Rscript
library(minfi)

# set the locations of important files
loc <- "/data/f114798/"
loc.idat <- "Data/ImageData"
loc.sheet <- "Samplesheets"
loc.rgdata <- "Data/RG_data"
loc.data <- "Data"
setwd(loc)
# conversion and storing of idat to rg files with option to subset on bad samples
targets <- read.metharray.sheet(loc.sheet, pattern="*.csv", recursive=TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
det.p <- detectionP(rg.set)
save(det.p, file=paste(loc.rgdata, "Detection_Pvalues.Rdata", sep="/"))
# remove sample if more then 20% of probes failed on the sample
bad.samples <- colMeans(det.p > 0.01) > 0.05
bad.sample.names.detP <- colnames(det.p[,bad.samples])
bad.sample.detp.rs <- data.frame(CallRate=1-colMeans(det.p[,bad.samples] > 0.01))
rownames(bad.sample.detp.rs) <- rg.set$Sample_Name[rg.set$Basename %in% bad.sample.names.detP]
write.csv(bad.sample.detp.rs, file=paste(loc.data, "Bad_Samples_CallRate.csv", sep="/"))
save(rg.set, targets, bad.sample.names.detP, bad.samples, file=paste(loc.rgdata, "RG_Channel_Set.Rdata", sep = "/"))
targets <- targets[!bad.samples,]
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(rg.set, targets, bad.sample.names.detP, bad.samples, file=paste(loc.rgdata, "RG_Channel_Set_Clean.Rdata", sep = "/"))

