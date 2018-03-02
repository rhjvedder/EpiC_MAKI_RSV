#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(gridExtra)

start.time <- Sys.time()

# run on peregrine or a system with atleast R 3.4.2 and minfi 3.6
# edit the output place for logs and output
sink("/data/f114798/Logs/log_job_p1.log")
loc.out <- "/data/f114798/Data"

# loading data
load(file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))
qc.total.probes <- nrow(rg.set)
qc.total.samples <- length(rg.set$Sample_Name)
manifest <- getManifest(rg.set)
print("manifest")
manifest
# location of phenotype data
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# Annotation
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# QC on signal detection removing samples with detection p-values over 0.05
print("QC on signal detection removing samples with detection p-values over 0.05")
det.p <- detectionP(rg.set)
save(det.p, file=paste(loc.out, "epic_maki_detP.Rdata", sep="/"))

# remove sample if more then 20% of probes failed on the sample
bad.samples <- colMeans(det.p > 0.01) > 0.05
bad.sample.names.detP <- colnames(det.p[,bad.samples])
rg.set <- rg.set[,!bad.samples]
targets <- targets[!bad.samples,]
all.det.p <- det.p
det.p <- det.p[,!bad.samples]
table(bad.samples)

# remove any probes that have failed in 50% of samples
bad.probes <- rowMeans(det.p > 0.01) > 0.1
table(bad.probes)
bad.probe.names.detP <- rownames(det.p[bad.probes,])

save(det.p, file=paste(loc.out, "epic_maki_detP.Rdata", sep="/"))
save(rg.set, targets, file=paste(loc.out, "epic_maki_rgSet_postDetP.Rdata", sep="/"))
#save 

# add sex to samples
m.set <- preprocessRaw(rg.set)
gm.set <- mapToGenome(m.set)
save(m.set, targets, file=paste(loc.out, "epic_maki_mSet.Rdata", sep="/"))
save(gm.set, targets, file=paste(loc.out, "epic_maki_gmSet.Rdata", sep="/"))

predicted.sex <- getSex(gm.set, cutoff = -2)
gm.set <- addSex(gm.set, predicted.sex$predictedSex)

# remove samples of which the reported sex does not equal the sex chromosomes
# keep only the phenotype data of dimthe samples
# prevent lowercase characters in data
gm.set$Sample_Name <- toupper(gm.set$Sample_Name)
samples.gm <- data.frame(Sample=gm.set$Sample_Name, Gender=gm.set$predictedSex, stringsAsFactors = FALSE)
samples.in <- data.frame(Sample=Phenotype[,2], Number=rep(1, length(Phenotype[,2])), Gender=Phenotype[,6], stringsAsFactors = FALSE)

save(gm.set, targets, file=paste(loc.out, "epic_maki_gmSet_wPS.Rdata", sep="/"))

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

# PC of X chromosomes
sex.chromosome.probes <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% "chrX"]
x.chromosome.probes <- gm.set[sex.chromosome.probes,]
x.beta <- getBeta(x.chromosome.probes)
x.beta <- na.omit(x.beta)
x.beta <- t(x.beta)
##  prop is % variance explained
samples.code <- sapply(rownames(x.beta), function(x) {as.character(gm.set$Sample_Name[gm.set$Basename==x])})
x.pheno <- sapply(samples.code, function(x) {
  gender <- samples.in$Gender[samples.in$Sample==x]
  if(!is.null(gender) && length(gender) == 1 && !is.na(gender)) {
    if (gender=="F") {
      y <- "F"
    } else if (gender=="M") {
      y <- "M"
    } else {
      y <- NA
    }
  } else {
    y <- NA
  }
  y
})
x.pca.Color <- sapply(x.pheno, function(x) {
  if (!is.na(x)) {
    if (x=="F") {
      y <- "red"
    } else if (x=="M") {
      y <- "blue"
    } else {
      y <- "grey"
    }
  } else {
    y <- "grey"
  }
  y
})
## the row of data should be samples. the default setting center=TRUE and scale=FALSE
out <- prcomp(x.beta)
##  change 1 or 2 , you can get different plots
prop<-out$sdev^2 / sum(out$sdev^2) 
png(paste(loc.out, "pca_gender_prcomp.png", sep = "/"))
plot(out$x[,1],out$x[,2],xlab="PC1",ylab="PC2", main="PCA of Gender", col=x.pca.Color)
print("---------------------------------------------------------------------")
dev.off()


bad.sample.genders <- x.pheno[rownames(x.beta) %in% bad.samples.selection]
bad.samples.data <- x.beta[rownames(x.beta) %in% bad.samples.selection,]
# bad.samples.data <- apply(x.beta, 1, bad.samples.data)

# Histogram of Y chromosomes to test posssible xxy
sex.chromosome.probes <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% "chrY"]
# print out which do not have y chromosome probes
y.chromosome.probes <- gm.set[sex.chromosome.probes,]
y.set <- y.chromosome.probes[,colnames(y.chromosome.probes) %in% bad.samples.selection]
y.beta <- getBeta(y.set)
y.beta <- na.omit(y.beta)
y.beta <- t(y.beta)

for (row in 1:nrow(bad.samples.data)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.samples.selection[row], "_", bad.sample.genders[row], "_histogram_x.png"), sep = "/"))
    hist(bad.samples.data[row,], main=paste("Histogram of ", bad.samples.selection[row], "'s X chromosome beta values"))
  dev.off()
}

for (row in 1:nrow(y.beta)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.samples.selection[row], "_", bad.sample.genders[row], "_histogram_y.png"), sep = "/"))
    hist(y.beta[row,], main=paste("Histogram of ", bad.samples.selection[row], "'s X chromosome beta values"))
  dev.off()
}

good.samples.names <- samples.gm$Sample[samples.good]
table(samples.good)
gm.set <- gm.set[,samples.good]
targets <- targets[samples.good,]
print("---------------------------------------------------------------------")

# remove probes with bad det p from gm set
gm.set <- gm.set[!bad.probes,]

# if your data includes males and females, remove probes on the sex chromosomes
print("if your data includes males and females, remove probes on the sex chromosomes")
keep <- !(featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
bad.probe.names.xy <- rownames(gm.set[!keep,])
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter cross reactive probes
print("filter cross reactive probes")
reactive.probes <- read.csv(file="/data/f114798/Rscripts/1-s2.0-S221359601630071X-mmc2.csv", sep="\t", stringsAsFactors=FALSE)
keep <- !(featureNames(gm.set) %in% reactive.probes$IlmnID)
bad.probe.names.cr <- rownames(gm.set[!keep,])
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")

# filter polymorphic targets
print("filter polymorphic probes")
polymorphic.probes <- read.csv(file="/data/f114798/Rscripts/1-s2.0-S221359601630071X-mmc1.csv", sep="\t", stringsAsFactors=FALSE)
## $EUR_AF < 0.05 remove
keep <- !(featureNames(gm.set) %in% polymorphic.probes$IlmnID[polymorphic.probes$EUR_AF<0.05])
bad.probe.names.pm <- rownames(gm.set[!keep,])
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")


# remove probes with SNPs at CpG site ; same as above
#print("remove probes with SNPs at CpG sites")
#gm.set.a <- dropLociWithSnps(gm.set, maf=0.05)
#keep <- rownames(gm.set) %in% rownames(gm.set.a)
#bad.probe.names.SC <- rownames(gm.set[!keep,])
#table(keep)
#gm.set <- gm.set.a
#print("---------------------------------------------------------------------")


# filtering step to take out certain SNPs starting with 'rs'
print("filtering step to take out certain SNPs starting with 'rs'")
keep <- !(grepl("rs", featureNames(gm.set)))
bad.probe.names.rs <- rownames(gm.set[!keep,])
gm.set <- gm.set[keep,]
table(keep)
print("---------------------------------------------------------------------")
# maybe change the order of filtering and normalizing

# original number of probes
print("original amount of probes")
qc.total.probes
# remaining number of probes
removed.probes <- qc.total.probes - nrow(gm.set)
remaining.probes <- qc.total.probes - removed.probes
print("remaining number of probes")
remaining.probes

# note original total amount of samples
qc.probes.colnames <- c("Detect p-value", "Probes X/Y", "Cross Reactive", "Polymorphic Pr", "Probes with 'rs'", "Probes left")
qc.samples.colnames <- c("Detect p-value", "Gender check", "Samples left")

nbr_sampl <- function(x) {return(if (is.null(x)) {0} else {length(x)})}
nam_sampl <- function(x) {return(if (is.null(x)) {"No Samples removed"} else {x})}

qc.s.df <- data.frame(amount_detP=nbr_sampl(bad.sample.names.detP), 
                      amount_Gender=nbr_sampl(bad.sample.names.gender),
                      amount_left=(qc.total.samples - nbr_sampl(bad.sample.names.detP)) - nbr_sampl(bad.sample.names.gender))
colnames(qc.s.df) <- qc.samples.colnames
pdf(paste(loc.out, "qc_Sample_Table.pdf", sep = "/"), height=11, width=8.5)
grid.table(qc.s.df)
dev.off()
print("---------------------------------------------------------------------")
print("removed on detection p-value")
print(nam_sampl(bad.sample.names.detP))
print("removed on gender")
print(nam_sampl(bad.sample.names.gender))
print("---------------------------------------------------------------------")
qc.p.df <- data.frame(nbr_sampl(bad.probe.names.detP), nbr_sampl(bad.probe.names.xy), nbr_sampl(bad.probe.names.cr), nbr_sampl(bad.probe.names.pm), nbr_sampl(bad.probe.names.rs), remaining.probes)
colnames(qc.p.df) <- qc.probes.colnames
pdf(paste(loc.out, "qc_Probe_Table.pdf", sep = "/"), height=11, width=8.5)
grid.table(qc.p.df)
dev.off()
bad.probe.names.all <- c(nam_sampl(bad.probe.names.detP), nam_sampl(bad.probe.names.xy), nam_sampl(bad.probe.names.cr), nam_sampl(bad.probe.names.pm), nam_sampl(bad.probe.names.rs))
sink()
sink("/data/f114798/Data/qc_bad_probes.csv")
print(bad.probe.names.all)
sink()
sink("/data/f114798/Logs/log_job_p2.log")

print("---------------------------------------------------------------------")
save(gm.set, targets, file = paste(loc.out, "epic_maki_filtered_gm_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
