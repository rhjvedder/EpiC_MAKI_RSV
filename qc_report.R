#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
loc.out <- "/data/f114798/Data"

load(file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))
load(file=paste(loc.out, "epic_maki_mSet.Rdata", sep="/"))
load(file=paste(loc.out, "epic_maki_gmSet.Rdata", sep="/"))

qc <- getQC(m.set)
png(paste(loc.out, "qc_plot.png", sep="/"))
plotQC(qc)
dev.off()

qcReport(rg.set, pdf= paste(loc.out, "qcReport.pdf", sep = "/"))
predicted.sex <- getSex(gm.set, cutoff = -2)
png(paste(loc.out, "predicted_sex.png", sep = "/"))
plotSex(predicted.sex)
dev.off()

