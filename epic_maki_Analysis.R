#!/usr/bin/env Rscript
library(minfi)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log_3.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Data"

# load data
load(file = paste(loc.out, "epic_maki_filtered_meth_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
print("manifest")
manifest
Phenotype <- read.csv("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# Analysis
beta <- getBeta(m.set.sq)
m <- getM(m.set.sq)

# remove outliers
remove_outliers <- function(x, na.rm = TRUE) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  extreme <- 3.0 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - extreme)] <- NA
  y[x > (qnt[2] + extreme)] <- NA
  y
}

M <- apply(m, 1, remove_outliers)
save(M, file = paste(loc.out, "M_Values_Maki.Rdata", sep = "/"))
# amount of extreme outliers
print("The amount of extreme outliers:")
which(is.na(M))
print("---------------------------------------------------------------------")

# idea for analyzing, show a histogram with the amount difference between probes

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
