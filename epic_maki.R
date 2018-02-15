#!/usr/bin/env Rscript
# installation of bioconductor, minfi and loading of libraries
# source("https://bioconductor.org/biocLite.R")
# install.packages("RColorBrewer")
# biocLite("minfi")
# biocLite("minfiDataEPIC")
library(minfi)
library(minfiData)
library(foreign)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# load in cml arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("This script needs an output dir, input dir and a normalizing method [ssNoob|Funnorm|Quantile].n", call.=FALSE)
}

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log.txt")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- args[1]
loc.data <- args[2]
norma = args[3]

# loading data
targets <- read.metharray.sheet(file.path(paste(loc.data, "SampleSheets", sep = "/")))
rg.set <- read.metharray.exp(file.path(paste(loc.data, "ImageData", sep = "/"), targets=targets))
save(rg.set, file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
print("manifest")
manifest
Phenotype <- read.csv("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/Age_at_MAKI_III_nasal_sampling_17_Jan_2018.csv")
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

# normalize based on the supplied normalizing type
if(norma == "Funnorm") {
  m.set.sq <- preprocessFunnorm(rg.set)
} else if(norma == "ssNoob") {
  m.set.sq <- preprocessNoob(rg.set, dyeMethod = "single")
  m.set.sq <- mapToGenome(m.set.sq)
} else {
  m.set.sq <- preprocessQuantile(rg.set)
}

# plots showing normalizing before and after
png(paste(loc.out, "normalizing.png", sep = "/"))
par(mfrow=c(1,2))
densityPlot(rg.set, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(m.set.sq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")

# ensure probes are in the same order in the m.set.sq and det.p objects
det.p <- det.p[match(featureNames(m.set.sq),rownames(det.p)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(det.p < 0.01) == ncol(m.set.sq)
m.set.sq.flt <- m.set.sq[keep,]
table(keep)

# add sex to samples
predicted.sex <- getSex(m.set.sq.flt, cutoff = -2)
png(paste(loc.out, "predicted_sex.png", sep = "/"))
plotSex(predicted.sex)
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")
m.set.sq.flt <- addSex(m.set.sq.flt, predicted.sex$predictedSex)

# remove samples of which the reported sex does not equal the sex chromosomes
# keep only the phenotype data of the samples
Phenotype <- Phenotype[Phenotype$Trailnummer %in% m.set.sq.flt$Sample_Name,]
# order the Phenotype data to the methyl set
Phenotype <- Phenotype[m.set.sq.flt$Sample_Name,]
# compare the phenotypes to the predicted data, which can indicate errors in lab
keep <- m.set.sq.flt$predictedSex == Phenotype$Geslacht
m.set.sq.flt <- m.set.sq.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# if your data includes males and females, remove probes on the sex chromosomes
print("if your data includes males and females, remove probes on the sex chromosomes")
keep <- !(featureNames(m.set.sq.flt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
m.set.sq.flt <- m.set.sq.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter cross reactive probes
print("filter cross reactive probes")
reactive.probes <- read.csv(file="/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/1-s2.0-S221359601630071X-mmc1.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(m.set.sq.flt) %in% reactive.probes$IlmnID)
m.set.sq.flt <- m.set.sq.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter polymorphic targets
print("filter polymorphic probes")
polymorphic.probes <- read.csv(file="/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/1-s2.0-S221359601630071X-mmc2.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(m.set.sq.flt) %in% polymorphic.probes$IlmnID)
m.set.sq.flt <- m.set.sq.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")


# remove probes with SNPs at CpG site
print("remove probes with SNPs at CpG sites")
m.set.sq.flt <- dropLociWithSnps(m.set.sq.flt, maf=0.05)
print("---------------------------------------------------------------------")


# filtering step to take out certain SNPs starting with 'rs'
print("filtering step to take out certain SNPs starting with 'rs'")
keep <- !(grepl("rs", featureNames(m.set.sq.flt)))
m.set.sq.flt <- m.set.sq.flt[keep,]
table(keep)
print("---------------------------------------------------------------------")
# maybe change the order of filtering and normalizing

# original number of probes
print("original amount of probes")
total.probes
# remaining number of probes
removed.probes <- total.probes - nrow(m.set.sq.flt)
remaining.probes <- total.probes - removed.probes
print("remaining number of probes")
remaining.probes
print("---------------------------------------------------------------------")
save(m.set.sq.flt, file = paste(loc.out, "epic_maki_filtered_meth.Rdata", sep = "/"))

# Analysis
beta <- getBeta(m.set.sq.flt)
m <- getM(m.set.sq.flt)

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
