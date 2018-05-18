#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
loc <- "/data/f114798/"
loc.rgdata <- "Data/RG_data"
loc.qc <- "Data/QC"
loc.data <- "Data/meth_set"
setwd(loc)

load(file = paste(loc.rgdata, "RG_Channel_Set.Rdata", sep = "/"))
load(file=paste(loc.data, "epic_maki_mSet.Rdata", sep="/"))
load(file=paste(loc.data, "epic_maki_gmSet.Rdata", sep="/"))

qc <- getQC(m.set)
png(paste(loc.qc, "qc_plot.png", sep="/"))
plotQC(qc)
dev.off()

qcReport(rg.set, pdf= paste(loc.qc, "qcReport.pdf", sep = "/"))
predicted.sex <- getSex(gm.set, cutoff = -2)
png(paste(loc.qc, "predicted_sex.png", sep = "/"))
plotSex(predicted.sex)
dev.off()

