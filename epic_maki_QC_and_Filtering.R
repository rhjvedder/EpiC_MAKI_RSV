#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gridExtra)

sink("/data/f114798/Logs/log_sex.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/data/f114798/Data"

# loading data
load(file = paste("/data/f114798/Data", "epic_maki_rg_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
print("manifesrt")
manifest
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# Annotation
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# note original total amount of samples
qc.probes.colnames <- c("Detect p-value", "Failed Probes", "Probes X/Y", "Cross Reactive", "Polymorphic Pr", "Probes SNPs@CpGs", "Probes with 'rs'")
qc.samples.colnames <- c("Detect p-value", "Gender check")
qc.total.probes <- nrow(rg.set)

# QC on signal detection removing samples with detection p-values over 0.05
print("QC on signal detection removing samples with detection p-values over 0.05")
det.p <- detectionP(rg.set)

# remove sample if more then 20% of probes failed on the sample
bad.samples <- colMeans(det.p > 0.01) > 0.05
removed.sample.names.detP <- colnames(det.p[,bad.samples])
rg.set <- rg.set[,!bad.samples]
targets <- targets[!bad.samples,]
det.p <- det.p[,!bad.samples]
table(bad.samples)

# remove any probes that have failed in 50% of samples
bad.probes <- rowMeans(det.p > 0.01) > 0.1
table(bad.probes)
bad.probes.names.detP <- rownames(det.p[bad.probes,])

save(det.p, file=paste(loc.out, "epic_maki_detP.Rdata", sep="/"))
save(rg.set, file=paste(loc.out, "epic_maki_rgSet_postDetP.Rdata", sep="/"))
#save 

# add sex to samples
m.set <- preprocessRaw(rg.set)
gm.set <- mapToGenome(m.set)
save(gm.set, file=paste(loc.out, "epic_maki_gmSet.Rdata", sep="/"))

qc <- getQC(m.set)
png(paste(loc.out, "qc_plot.png", sep="/"))
plotQC(qc)
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")
qcReport(rg.set, pdf= paste(loc.out, "qcReport.pdf", sep = "/"))

predicted.sex <- getSex(gm.set, cutoff = -2)
png(paste(loc.out, "predicted_sex.png", sep = "/"))
plotSex(predicted.sex)
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")
gm.set <- addSex(gm.set, predicted.sex$predictedSex)

# remove samples of which the reported sex does not equal the sex chromosomes
# keep only the phenotype data of dimthe samples
# prevent lowercase characters in data
gm.set$Sample_Name <- toupper(gm.set$Sample_Name)
samples.gm <- data.frame(Sample=gm.set$Sample_Name, Gender=gm.set$predictedSex, stringsAsFactors = FALSE)
samples.in <- data.frame(Sample=Phenotype[,2], Number=rep(1, length(Phenotype[,2])), Gender=Phenotype[,6], stringsAsFactors = FALSE)

samples.good <- apply(samples.gm, 1, function(x) {
  if (x[1] %in% samples.in$Sample) {
    i <- samples.in$Number[samples.in$Sample==x[1]]
    y <- samples.in$Gender[samples.in$Sample==x[1]]
    if (i > 1) {
      z <- x[2] == y[i]
    } else {
      z <- x[2] == y
    }
    samples.in[samples.in$Sample==x[1],2] <- samples.in$Number[samples.in$Sample==x[1]] + 1
  } else {
    z <- TRUE
  }
  z
})

bad.samples.names <- samples.gm$Sample[!samples.good]
good.samples.names <- samples.gm$Sample[samples.good]
table(samples.good)
gm.set <- gm.set[,samples.good]
targets <- targets[samples.good,]
print("---------------------------------------------------------------------")

# PC of X chromosomes
sex.chromosome.probes <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% "chrX")]
x.chromosome.probes <- gm.set[sex.chromosome.probes,]
x.beta <- getBeta(x.chromosome.probes)
x.beta <- na.omit(x.beta)
x.beta <- t(x.beta)
x.covmat <- cov(x.beta)
x.eig <- eigen(x.covmat)



png(paste(loc.out, "predicted_sex_postGenderCheck.png", sep = "/"))
plotSex(gm.set$predictedSex)
print("---------------------------------------------------------------------")
dev.off()

# if your data includes males and females, remove probes on the sex chromosomes
print("if your data includes males and females, remove probes on the sex chromosomes")
keep <- !(featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter cross reactive probes
print("filter cross reactive probes")
reactive.probes <- read.csv(file="/data/f114798/Rscripts/1-s2.0-S221359601630071X-mmc1.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(gm.set) %in% reactive.probes$IlmnID)
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter polymorphic targets
print("filter polymorphic probes")
polymorphic.probes <- read.csv(file="/data/f114798/Rscripts/1-s2.0-S221359601630071X-mmc2.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(gm.set) %in% polymorphic.probes$IlmnID)
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")


# remove probes with SNPs at CpG site
print("remove probes with SNPs at CpG sites")
gm.set.a <- dropLociWithSnps(gm.set, maf=0.05)
keep <- gm.set %in% gm.set.a
table(keep)
gm.set < gm.set.a
print("---------------------------------------------------------------------")


# filtering step to take out certain SNPs starting with 'rs'
print("filtering step to take out certain SNPs starting with 'rs'")
keep <- !(grepl("rs", featureNames(gm.set)))
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")
# maybe change the order of filtering and normalizing

# original number of probes
print("original amount of probes")
total.probes
# remaining number of probes
removed.probes <- total.probes - nrow(gm.set)
remaining.probes <- total.probes - removed.probes
print("remaining number of probes")
remaining.probes


pdf("qcTable.pdf", height=11, width=8.5)
grid.table(qc.df)
dev.off()

print("---------------------------------------------------------------------")
rg.set.flt <- rg.set.flt[rg.set.flt %in% gm.set]
save(rg.set.flt, targets, file = paste(loc.out, "epic_maki_filtered_rg_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
