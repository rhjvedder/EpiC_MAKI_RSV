#!/usr/bin/env Rscript
source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylationEPICmanifest")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log_1.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Data"

# loading data
load(file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
print("manifest")
manifest
Phenotype <- read.csv("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# Annotation
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# note original total amount of samples
total.probes <- nrow(rg.set)

# QC on signal detection removing samples with detection p-values over 0.05
print("QC on signal detection removing samples with detection p-values over 0.05")
det.p <- detectionP(rg.set)
keep <- colMeans(det.p) < 0.05
rg.set <- rg.set[,keep]
targets <- targets[keep,]
det.p <- det.p[,keep]
table(keep)

# ensure probes are in the same order in the m.set.sq and det.p objects
det.p <- det.p[match(featureNames(rg.set),rownames(det.p)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(det.p < 0.01) == ncol(rg.set)
rg.set.flt <- rg.set[keep,]
table(keep)

# add sex to samples
predicted.sex <- getSex(rg.set.flt, cutoff = -2)
png(paste(loc.out, "predicted_sex.png", sep = "/"))
plotSex(predicted.sex)
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")
rg.set.flt <- addSex(rg.set.flt, predicted.sex$predictedSex)

# remove samples of which the reported sex does not equal the sex chromosomes
# keep only the phenotype data of the samples
Phenotype.avail <- Phenotype[Phenotype$Trailnummer %in% rg.set.flt$Sample_Name,]
# order the Phenotype data to the methyl set
Phenotype.avail <- Phenotype.avail[rg.set.flt$Sample_Name,]
# compare the phenotypes to the predicted data, which can indicate errors in lab
keep <- rg.set.flt$predictedSex == Phenotype.avail$Geslacht
rg.set.flt <- rg.set.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# if your data includes males and females, remove probes on the sex chromosomes
print("if your data includes males and females, remove probes on the sex chromosomes")
keep <- !(featureNames(rg.set.flt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
rg.set.flt <- rg.set.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter cross reactive probes
print("filter cross reactive probes")
reactive.probes <- read.csv(file="/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/1-s2.0-S221359601630071X-mmc1.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(rg.set.flt) %in% reactive.probes$IlmnID)
rg.set.flt <- rg.set.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter polymorphic targets
print("filter polymorphic probes")
polymorphic.probes <- read.csv(file="/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/1-s2.0-S221359601630071X-mmc2.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(rg.set.flt) %in% polymorphic.probes$IlmnID)
rg.set.flt <- rg.set.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")


# remove probes with SNPs at CpG site
print("remove probes with SNPs at CpG sites")
rg.set.flt <- dropLociWithSnps(rg.set.flt, maf=0.05)
print("---------------------------------------------------------------------")


# filtering step to take out certain SNPs starting with 'rs'
print("filtering step to take out certain SNPs starting with 'rs'")
keep <- !(grepl("rs", featureNames(rg.set.flt)))
rg.set.flt <- rg.set.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")
# maybe change the order of filtering and normalizing

# original number of probes
print("original amount of probes")
total.probes
# remaining number of probes
removed.probes <- total.probes - nrow(rg.set.flt)
remaining.probes <- total.probes - removed.probes
print("remaining number of probes")
remaining.probes
print("---------------------------------------------------------------------")
save(rg.set.flt, file = paste(loc.out, "epic_maki_filtered_rg_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
