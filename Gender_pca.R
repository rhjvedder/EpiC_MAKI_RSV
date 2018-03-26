#!/usr/bin/env Rscript
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# example data EpiC 850K
loc.out <- "/data/f114798/Data"
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
load(file=paste(loc.out, "epic_maki_gmSet_wPS.Rdata", sep="/"))
load(file=paste(loc.out, "epic_maki_gmSet_wPS.Rdata", sep="/"))

# PC of X chromosomes
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

sex.chromosome.probes <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% "chrX"]
x.chromosome.probes <- gm.set[sex.chromosome.probes,]
x.beta <- getBeta(x.chromosome.probes)
x.beta <- na.omit(x.beta)
x.beta <- apply(x.beta, 1, scale)
x.beta <- t(x.beta)
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
# use pca from chen and compare to a scaled beta
x.covmat <- cov(x.beta)
x.eig <- eigen(x.covmat)
100*cumsum(x.eig$values[1:15])/sum(x.eig$values)
prj <- as.matrix(x.beta) %*% x.eig$vectors
x.pca <- data.frame(PC1=prj[,1], PC2=prj[,2], PG=x.pheno)
x.pca$Color <- sapply(x.pheno, function(x) {
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
png(paste(loc.out, "pca_gender.png", sep = "/"))
plot(x.pca$PC1,x.pca$PC2,xlab="PC1", ylab="PC2", main="PCA of Gender", col=x.pca$Color)
print("---------------------------------------------------------------------")
dev.off()

bad.samples.names.gender <- samples.gm$Sample[!samples.good]
print(bad.samples.names.gender)

bad.samples.selection <- gm.set$Basename[!samples.good]
bad.samples.data <- x.beta[rownames(x.beta) %in% bad.samples.selection,]
# bad.samples.data <- apply(x.beta, 1, bad.samples.data)

for (row in 1:nrow(bad.samples.data)) {
  png(paste(loc.out, "Gender_histograms", paste(bad.samples.selection[row], "_histogram.png"), sep = "/"))
    hist(bad.samples.data[row,], main=paste("Histogram of ", bad.samples.selection[row], "'s X chromosome beta values"))
  dev.off()
}
