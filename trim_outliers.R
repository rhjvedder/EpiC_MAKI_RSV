#!/usr/bin/env Rscript
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/data/f114798/Logs/log_trim_outliers.log")

# example data EpiC 850K
loc.out <- "/data/f114798/Data"

# load data
load(file=paste(loc.out, "epic_maki_ Quantile _meth_set.Rdata", sep = "/"))
load(file=paste(loc.out, "epic_maki_filtered_gm_set.Rdata", sep = "/"))
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
targets$Sample_Plate <- trimws(targets$Sample_Plate, "r")
targets$Sample_Name <- as.vector(sapply(targets$Sample_Name, function(x) {
  if (x!="Control_M") {
    y <- substr(x, 1, 4)
  } else {
    y <- x
  }
  y <- toupper(y)
  y
}, simplify=T))

m.set.sq$Sample_Name <- as.vector(sapply(m.set.sq$Sample_Name, function(x) {
  if (x!="Control_M") {
    y <- substr(x, 1, 4)
  } else {
    y <- x
  }
  y <- toupper(y)
  y
}, simplify=T))

m.set.flt <- m.set.sq[rownames(m.set.sq) %in% rownames(gm.set),colnames(m.set.sq) %in% colnames(gm.set)]
m.set.flt <- m.set.flt[,m.set.flt$Sample_Name!="CONTROL_M"]
targets <- targets[targets$Sample_Name %in% m.set.flt$Sample_Name,]

save(m.set.flt, file=paste(loc.out, "Methylation_Set.Rdata", sep="/"))
setwd(loc.out)
M.val <- getM(m.set.flt)
samplenames <- colnames(m.set.flt)
removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

system.time(OutlierResults<-removeOutliers(M.val))
M.val2<-OutlierResults[[1]]
Log<-OutlierResults[[2]]
save(M.val, targets, file="MAKI_trimmed_M.Rdata")
save(Log,file="Outlier_log_M.Rdata")

sink()
