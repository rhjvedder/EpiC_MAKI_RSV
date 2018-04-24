#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# This file contains functions for the pipeline to process 850k epic illumina data
# loading the rg.set
loc.idat <- "/data/f114798/Data/ImageData"
loc.sheet <- "/data/f114798/Samplesheets"
loc.out <- "/data/f114798/Data"

start.time <- Sys.time()

sink("/data/f114798/Logs/log_qc_p1.log")

# loading data
load(file=paste(loc.out, "Detection_Pvalues.Rdata", sep="/"))
load(file = paste(loc.out, "RG_Channel_Set.Rdata", sep = "/"))
qc.total.samples <- length(rg.set$Sample_Name)
manifest <- getManifest(rg.set)
print("manifest")
manifest
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# Annotation
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


# QC on signal detection removing samples with detection p-values over 0.05
print("QC on signal detection removing samples with detection p-values over 0.05")
all.det.p <- det.p
det.p <- det.p[,bad.samples]

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
qc.total.probes <- nrow(gm.set)
length(rownames(gm.set))
length(unique(rownames(gm.set)))
save(m.set, targets, file=paste(loc.out, "epic_maki_mSet.Rdata", sep="/"))
save(gm.set, targets, file=paste(loc.out, "epic_maki_gmSet.Rdata", sep="/"))

predicted.sex <- getSex(gm.set, cutoff = -2)
gm.set <- addSex(gm.set, predicted.sex$predictedSex)

# remove samples of which the reported sex does not equal the sex chromosomes
# keep only the phenotype data of dimthe samples
# prevent lowercase characters in data
gm.set$Sample_Name <- toupper(gm.set$Sample_Name)
samples.gm <- data.frame(Sample=gm.set$Sample_Name, Gender=gm.set$predictedSex, Basename=gm.set$Basename, stringsAsFactors = FALSE)
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
print("amount of X chromosome probes")
sum(sex.chromosome.probes == TRUE)
x.chromosome.probes <- gm.set[sex.chromosome.probes,]
x.beta <- getBeta(x.chromosome.probes)
x.beta <- na.omit(x.beta)
x.beta <- t(x.beta)
##  prop is % variance explained
samples.code <- sapply(rownames(x.beta), function(x) {as.character(gm.set$Sample_Name[gm.set$Basename==x])})
# add a color based on the gender, NAs are to be replaced with blue as these are controls and arent in the phenotype data
x.pheno <- sapply(samples.code, function(x) {
  gender <- samples.in$Gender[samples.in$Sample==x]
  if(!is.null(gender) && length(gender) == 1 && !is.na(gender)) {
    if (gender=="F") {
      y <- "F"
    } else if (gender=="M") {
      y <- "M"
    } else {
      y <- "M"
    }
  } else {
    y <- "M"
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
pca.out <- prcomp(x.beta)
##  change 1 or 2 , you can get different plots
prop<-out$sdev^2 / sum(out$sdev^2) 
png(paste(loc.out, "pca_x_gender_prcomp.png", sep = "/"))
ggplot(pca.out$x,aes(x=PC1,y=PC2,col=Phenotype.samples$Gender)) + geom_point(size=3,alpha=0.5) + labs(title="PCA plot of Gender, no controls, X chromosomes", x=paste("PC1 (", round(prop[1]*100, 2), "%)"), y=paste("PC2 (", round(prop[2]*100, 2), "%)"), colour="Gender")
dev.off()

# select out any samples which aren't grouped in with male or female
bad.quality.selection <- rownames(x.beta)[(out$x[,1] > -4) & (out$x[,1] < 10)]
bad.quality.trail <- gm.set$Sample_Name[gm.set$Basename %in% bad.quality.selection]
# bad quality samples?
bad.quality.trail
# make histograms for the four which are not misgendered

# make histograms of bad X chromosome methylation patterns of b-values
bad.samples.selection <- gm.set$Basename[!samples.good]
bad.sample.genders <- x.pheno[rownames(x.beta) %in% bad.samples.selection]
bad.samples.data <- x.beta[rownames(x.beta) %in% bad.samples.selection,]
bad.samples.trail <- gm.set$Sample_Name[!samples.good]
# bad.samples.data <- apply(x.beta, 1, bad.samples.data)

# Histogram of Y chromosomes to test posssible xxy
sex.chromosome.probes <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% "chrY"]
print("amount of Y chromosome probes")
sum(sex.chromosome.probes == TRUE)
# print out which do not have y chromosome probes
y.chromosome.probes <- gm.set[sex.chromosome.probes,]
y.set <- y.chromosome.probes[,colnames(y.chromosome.probes) %in% bad.samples.selection]
y.beta <- getBeta(y.set)
y.beta <- na.omit(y.beta)
y.beta <- t(y.beta)

# print x chromosome
for (row in 1:nrow(bad.samples.data)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.samples.trail[row], "_", bad.sample.genders[row], "_histogram_x.png"), sep = "/"))
    hist(bad.samples.data[row,], main=paste("Histogram of ", bad.sample.genders[row]," ", bad.samples.trail[row], "'s X chromosome beta values"), xlab="probe beta values")
  dev.off()
}

# print y chromosome
for (row in 1:nrow(y.beta)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.samples.trail[row], "_", bad.sample.genders[row], "_histogram_y.png"), sep = "/"))
    hist(y.beta[row,], main=paste("Histogram of ", bad.sample.genders[row], " ", bad.samples.trail[row], "'s Y chromosome beta values"), xlab="probe beta values")
  dev.off()
}

# print out the x and y chromosome histograms for methylation patterns by b-values
good.samples.names <- samples.gm$Basename[samples.good]
good.samples.genders <- x.pheno[rownames(x.beta) %in% good.samples.names]
good.samples.selection <- good.samples.names[min(which(good.samples.genders == "M"))]
good.samples.selection <- c(good.samples.names[min(which(good.samples.genders == "F"))], good.samples.selection)
good.samples.genders <- x.pheno[rownames(x.beta) %in% good.samples.selection]
good.samples.data <- x.beta[rownames(x.beta) %in% good.samples.selection,]
good.samples.trail <- gm.set$Sample_Name[gm.set$Basename %in% good.samples.selection]

# example x and y for good samples
y.set <- y.chromosome.probes[,colnames(y.chromosome.probes) %in% good.samples.selection]
y.beta <- getBeta(y.set)
y.beta <- na.omit(y.beta)
y.beta <- t(y.beta)

for (row in 1:nrow(good.samples.data)) {
  png(paste(loc.out, "Gender_histograms", paste("G_", good.samples.trail[row], "_", good.samples.genders[row], "_histogram_x.png"), sep = "/"))
    hist(good.samples.data[row,], main=paste("Histogram of good sample ", good.samples.genders[row]," ", good.samples.trail[row], "'s X chromosome beta values"), xlab="probe beta values")
  dev.off()
}

for (row in 1:nrow(y.beta)) {
  png(paste(loc.out, "Gender_histograms", paste("G_", good.samples.trail[row], "_", good.samples.genders[row], "_histogram_y.png"), sep = "/"))
    hist(y.beta[row,], main=paste("Histogram of good sample ", good.samples.genders[row], " ", good.samples.trail[row], "'s Y chromosome beta values"), xlab="probe beta values")
  dev.off()
}

bad.quality.trail.non_mg <- bad.quality.trail[bad.quality.trail %in% bad.samples.trail]
bad.quality.genders <- x.pheno[rownames(x.beta) %in% bad.quality.selection]
bad.quality.data <- x.beta[rownames(x.beta) %in% bad.quality.selection,]
# bad.samples.data <- apply(x.beta, 1, bad.samples.data)

y.set <- y.chromosome.probes[,colnames(y.chromosome.probes) %in% bad.quality.selection]
y.beta <- getBeta(y.set)
y.beta <- na.omit(y.beta)
y.beta <- t(y.beta)

# print x chromosome
for (row in 1:nrow(bad.quality.data)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.quality.trail[row], "_", bad.sample.genders[row], "_histogram_x.png"), sep = "/"))
  hist(bad.quality.data[row,], main=paste("Histogram of ", bad.sample.genders[row]," ", bad.quality.trail[row], "'s X chromosome beta values"), xlab="probe beta values")
  dev.off()
}

# print y chromosome
for (row in 1:nrow(y.beta)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.quality.trail[row], "_", bad.sample.genders[row], "_histogram_y.png"), sep = "/"))
  hist(y.beta[row,], main=paste("Histogram of ", bad.sample.genders[row], " ", bad.quality.trail[row], "'s Y chromosome beta values"), xlab="probe beta values")
  dev.off()
}

# graphing of the mean read scores for the samples

# edit to get read scores for both bad pca center and other six bad quality samples
bad.quality.detp <- all.det.p[,which(colnames(all.det.p) %in% bad.quality.selection)]
bad.quality.score <- apply(bad.quality.detp, 2, function(x) {mean(x > 0.01)})
bad.quality.rs <- data.frame(CallRate=bad.quality.score)
rownames(bad.quality.rs) <- bad.quality.trail
write.csv(bad.quality.rs, file=paste(loc.out, "pca_center_RS.csv", sep="/"))

bad.samples.phenotypes <- Phenotype[Phenotype[,2] %in% bad.samples.trail,]
bad.samples.detp <- all.det.p[,which(colnames(all.det.p) %in% bad.samples.selection)]
bad.samples.score <- apply(bad.samples.detp, 2, function(x) {mean(x > 0.01)})
bad.samples.rs <- data.frame(CallRate=bad.samples.score)
rownames(bad.samples.rs) <- bad.samples.trail
bad.samples.rs <- data.frame(CallRate=apply(bad.samples.phenotypes, 1, function(x) {
  bad.samples.rs[x[2],]
}))
rownames(bad.samples.rs) <- bad.samples.phenotypes[,2]
bad.samples.asthma <- data.frame(Samples=bad.samples.phenotypes[,2], Gender=bad.samples.phenotypes[,6], Asthma=bad.samples.phenotypes[,20], Diagnosis=bad.samples.phenotypes[,21], CallRate=bad.samples.rs[,1])
write.csv(bad.samples.asthma, file=paste(loc.out, "bad_samples_asthma.csv", sep="/"))

png(paste(loc.out, "bad_Samples_asthma_Table.png", sep = "/"), width = 800, height = 400, units = "px")
grid.table(bad.samples.asthma)
dev.off()

table(samples.good)
gm.set <- gm.set[,samples.good]
targets <- targets[samples.good,]
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(rg.set, targets, file=paste(loc.out, "RG_Channel_Set.Rdata", sep = "/"))
m.set <- preprocessRaw(rg.set)
gm.set <- mapToGenome(m.set)
print("---------------------------------------------------------------------")


# remove probes with bad det p from gm set
samples <- gm.set$Basename[samples.good]
# maybe change the order of filtering and normalizing

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
keep <- !(featureNames(gm.set) %in% polymorphic.probes$IlmnID[polymorphic.probes$EUR_AF>0.05])
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

# original number of probes
print("original amount of probes")
qc.total.probes
# remaining number of probes
removed.probes <- qc.total.probes - nrow(gm.set)
remaining.probes <- nrow(gm.set)
print("remaining number of probes")
remaining.probes

# note original total amount of samples
qc.probes.colnames <- c("Detect p-value", "Probes X/Y", "Cross Reactive", "Polymorphic Pr", "Probes left")
qc.samples.colnames <- c("Detect p-value", "Gender check", "Samples left")

nbr_sampl <- function(x) {return(if (is.null(x)) {0} else {length(unique(x))})}
nam_sampl <- function(x) {return(if (is.null(x)) {"No Samples removed"} else {x})}

qc.s.df <- data.frame(amount_detP=nbr_sampl(bad.sample.names.detP), 
                      amount_Gender=nbr_sampl(bad.sample.genders),
                      amount_left=(qc.total.samples - nbr_sampl(bad.sample.names.detP)) - nbr_sampl(bad.sample.genders))
colnames(qc.s.df) <- qc.samples.colnames
pdf(paste(loc.out, "qc_Sample_Table.pdf", sep = "/"), height=11, width=8.5)
grid.table(qc.s.df)
dev.off()
print("---------------------------------------------------------------------")
print("removed on detection p-value")
print(nam_sampl(bad.sample.names.detP))
print("removed on gender")
print(nam_sampl(bad.sample.genders))
print("---------------------------------------------------------------------")
qc.p.df <- data.frame(nbr_sampl(bad.probe.names.detP), nbr_sampl(bad.probe.names.xy), nbr_sampl(bad.probe.names.cr), nbr_sampl(bad.probe.names.pm), remaining.probes)
qc.n.df <- list(DetectP=bad.probe.names.detP, X=rownames(x.chromosome.probes), Y=rownames(y.chromosome.probes), CrossReact=bad.probe.names.cr, PolyMorp=bad.probe.names.pm)
for (d in names(qc.n.df)) {
  n <- data.frame(qc.n.df[d])
  write.csv(n, file=paste(loc.out, paste(d, "filtered_probes.csv", sep="_"), sep="/"))
}
colnames(qc.p.df) <- qc.probes.colnames
pdf(paste(loc.out, "qc_Probe_Table.pdf", sep = "/"), height=11, width=8.5)
grid.table(qc.p.df)
dev.off()
bad.probe.names.all <- c(nam_sampl(bad.probe.names.detP), nam_sampl(bad.probe.names.xy), nam_sampl(bad.probe.names.cr), nam_sampl(bad.probe.names.pm))
sink()
sink("/data/f114798/Data/qc_bad_probes.csv")
print(bad.probe.names.all)
sink()
sink("/data/f114798/Logs/log_qc_p2.log")

print("---------------------------------------------------------------------")
save(gm.set, targets, file = paste(loc.out, "epic_maki_filtered_gm_set.Rdata", sep = "/"))

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()