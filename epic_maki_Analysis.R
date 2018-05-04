#!/usr/bin/env Rscript
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log_3.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/data/f114798/Data"

# load data
load(file=paste(loc.out, "epic_maki_Quantile_meth_set.Rdata", sep = "/"))
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

Phenotype.samples <- apply(targets, 1, function(x) {
  if (x[1]!="CONTROL_M") {
    # get Phenotype data and batch data
    y <- c(x[1], x[8], x[3], if(Phenotype[Phenotype[,2]==x[1],6]=="F") {"F"} else {"M"}, Phenotype[Phenotype[,2]==x[1],20]=="Wheeze <12 maanden of astmamedicatie gebruik", Phenotype[Phenotype[,2]==x[1],21]=="Asthma or BHR diagnosis physician")
  } else {
    y <- c(x[1], x[8], x[3], "M", "NA", "NA")
  }
  y
})

m.set <- m.set.sq[,m.set.sq$Sample_Name!="CONTROL_M"]
m.set.M <- getM(m.set)
m.set.M <- na.omit(m.set.M)
m.set.M <- t(m.set.M)
pca.out <- prcomp(m.set.M)
prop<-pca.out$sdev^2 / sum(pca.out$sdev^2)
Phenotype.samples <- Phenotype.samples[,Phenotype.samples[1,]!="CONTROL_M"]
Phenotype.samples <- data.frame(Sample=Phenotype.samples[1,], Basename=Phenotype.samples[2,], Plate=Phenotype.samples[3,], Gender=Phenotype.samples[4,], Wheeze=Phenotype.samples[5,], Asthma=Phenotype.samples[6,])

png(paste(loc.out, "pca_all_gender_prcomp.png", sep = "/"))
ggplot(pca.out$x,aes(x=PC1,y=PC2,col=Phenotype.samples$Gender)) + geom_point(size=3,alpha=0.5) + labs(title="PCA plot of Gender, no controls", x=paste("PC1 (", round(prop[1]*100, 2), "%)"), y=paste("PC2 (", round(prop[2]*100, 2), "%)"), colour="Gender")
dev.off()

png(paste(loc.out, "pca_all_plate_prcomp.png", sep = "/"))
ggplot(pca.out$x,aes(x=PC1,y=PC2,col=Phenotype.samples$Plate)) + geom_point(size=3,alpha=0.5) + labs(title="PCA plot of Plate, no controls", x=paste("PC1 (", round(prop[1]*100, 2), "%)"), y=paste("PC2 (", round(prop[2]*100, 2), "%)"), colour="Plate")
dev.off()

png(paste(loc.out, "pca_all_wheeze_prcomp.png", sep = "/"))
ggplot(pca.out$x,aes(x=PC1,y=PC2,col=Phenotype.samples$Wheeze)) + geom_point(size=3,alpha=0.5) + labs(title="PCA plot of Wheeze, no controls", x=paste("PC1 (", round(prop[1]*100, 2), "%)"), y=paste("PC2 (", round(prop[2]*100, 2), "%)"), colour="Wheeze")
dev.off()

png(paste(loc.out, "pca_all_asthma_prcomp.png", sep = "/"))
ggplot(pca.out$x,aes(x=PC1,y=PC2,col=Phenotype.samples$Asthma)) + geom_point(size=3,alpha=0.5) + labs(title="PCA plot of Asthma, no controls", x=paste("PC1 (", round(prop[1]*100, 2), "%)"), y=paste("PC2 (", round(prop[2]*100, 2), "%)"), colour="Asthma")
dev.off()
# remove control samples before pca!

# remove the bad probes
m.set.flt <- m.set.sq[rownames(m.set.sq) %in% rownames(gm.set),colnames(m.set.sq) %in% colnames(gm.set)]
m.set.flt <- m.set.flt[,m.set.flt$Sample_Name!="CONTROL_M"]
print("---------------------------------------------------------------------")
# remove control samples
# outlier removal

Sys.time() - start.time
# stop logging
sink()
