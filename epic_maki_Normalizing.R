#!/usr/bin/env Rscript
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/data/f114798/Logs/log_norm.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/data/f114798/Data"
loc.out.norm <- paste(loc.out, "Norm", sep="/")

# loading data
load(file = paste(loc.out, "RG_Channel_Set.Rdata", sep = "/"))
print("---------------------------------------------------------------------")
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
  save(m.set.sq, file = paste(loc.out, paste("epic_maki_", norma, "_meth_set.Rdata"), sep = "/"))
  return(m.set.sq)
})

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
