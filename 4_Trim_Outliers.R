#!/usr/bin/env Rscript
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# set the locations of important files
loc <- "/data/f114798/"
loc.rgdata <- "Data/RG_data"
loc.qc <- "Data/QC"
loc.comp <- "Complementary_data"
loc.data <- "Data/meth_set"
loc.mdata <- "Data/M_Values"
setwd(loc)

# load data
load(file=paste(loc.data, "epic_maki_Quantile_meth_set.Rdata", sep = "/"))
load(file=paste(loc.rgdata, "RG_Channel_Set_Clean.Rdata", sep = "/"))
probes <- read.csv(paste(loc.comp, 'filtered_probes.csv', sep = "/"), stringsAsFactors = FALSE)
probes <- c(probes[,1], probes[,2], probes[,3], probes[,4], probes[,5])
probes <- unique(probes)
probes <- probes[!(probes=="")]
l.p <- length(probes)
l.m <- length(rownames(rg.set))
l.m - l.p

Phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
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

m.set.flt <- m.set.sq[!(rownames(m.set.sq) %in% probes),]
print("before control removal")
dim(m.set.flt)
m.set.flt <- m.set.flt[,m.set.flt$Sample_Name!="CONTROL_M"]
print("after control removal")
dim(m.set.flt)
print("remaining")
length(rownames(m.set.flt))
targets <- targets[targets$Sample_Name %in% m.set.flt$Sample_Name,]
dim(m.set.flt)

save(m.set.flt, file=paste(loc.data, "Methylation_Set.Rdata", sep="/"))
M.val <- getM(m.set.flt)
B.val <- getBeta(m.set.flt)
print("amount of duplicated data")
colnames(M.val[,duplicated(targets$Sample_Name)])
M.val <- M.val[,!(duplicated(targets$Sample_Name))]
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
M.val<-OutlierResults[[1]]
print("data size after trimming")
dim(M.val)
Log<-OutlierResults[[2]]
save(M.val, targets, file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))
save(B.val, targets, file=paste(loc.mdata, "MAKI_B_vals.Rdata", sep="/"))
save(Log, file=paste(loc.mdata, "Outlier_log_M.Rdata", sep="/"))
