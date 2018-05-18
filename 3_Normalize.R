#!/usr/bin/env Rscript
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# set the locations of important files
loc <- "/data/f114798/"
loc.idat <- "Data/ImageData"
loc.sheet <- "Samplesheets"
loc.rgdata <- "Data/RG_data"
loc.data <- "Data/meth_set"
setwd(loc)

# load the samples to remove before normalization
samples <- read.csv(paste(loc.qc, 'bad_samples.csv', sep="/"), stringsAsFactors = FALSE)
samples <- unique(samples[,2])
# remove bad samples
targets <- read.metharray.sheet(loc.sheet, pattern="*.csv", recursive=TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})

targets <- targets[!(targets$Basename %in% samples),]
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(rg.set, targets, samples, file=paste(loc.rgdata, "RG_Channel_Set_Clean.Rdata", sep = "/"))

targets$Sample_Plate <- trimws(targets$Sample_Plate, "r")
# do a selection on the original rg set with the samples and probes from the gm set
# requires a selection on the original probe names
# rg.set.flt <- rg.set[rownames(rg.set) %in% rownames(gm.set),colnames(rg.set) %in% colnames(gm.set)]

norms <- c("Funnorm", "ssNoob", "Quantile")

norms.sets <- sapply(norms, function(norma) {
  # normalize based on the supplied normalizing type
  if(norma == "Funnorm") {
    m.set.sq <- preprocessFunnorm(rg.set)
  } else if(norma == "ssNoob") {
    m.set.sq <- preprocessNoob(rg.set, dyeMethod = "single")
    m.set.sq <- mapToGenome(m.set.sq)
  } else {
    m.set.sq <- preprocessQuantile(rg.set)
  }
  save(m.set.sq, file = paste(loc.data, paste0("epic_maki_", norma, "_meth_set.Rdata"), sep = "/"))
  return(m.set.sq)
})
