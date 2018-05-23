#!/usr/bin/env Rscript
library(WGCNA)
library(plyr)
library(parallel)
library(matrixStats)
options(stringsAsFactors = FALSE)
# set the locations of important files
loc <- "/data/f114798/"
loc.comp <- "Complementary_data"
loc.mdata <- "Data/M_Values"
loc.wcna <- "Data/Wcna"
setwd(loc)
load(file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))

phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
phenotype <- data.frame(Sample=phenotype[,2], Wheeze=phenotype[,20], Asthma=phenotype[,21], Prevention=phenotype[,27], FEV05=phenotype[,26], Gender=phenotype[,6], Age=phenotype[,28], IgE=phenotype[,19], Prednison=phenotype[,22], SmokingBirthMom=phenotype[,9], SmokingPregnacy=phenotype[,10], SmokingInsideBirth=phenotype[,12], stringsAsFactors=FALSE)
n <- length(targets$Sample_Name)
phenotype1 <- data.frame(Basename=rep(NA, n), Sample=rep(NA, n), Wheeze=rep(NA, n), Asthma=rep(NA, n), Prevention=rep(NA, n), FEV05=rep(NA, n), Gender=rep(NA, n), Age=rep(NA, n), IgE=rep(NA, n), Prednison=rep(NA, n), SmokingBirthMom=rep(NA, n), SmokingPregnacy=rep(NA, n), SmokingInsideBirth=rep(NA, n), stringsAsFactors=FALSE)

phenotype$Wheeze <- as.factor(phenotype$Wheeze)
phenotype$Asthma <- as.factor(phenotype$Asthma)
phenotype$Prevention <- as.factor(phenotype$Prevention)
phenotype$Gender <- as.factor(phenotype$Gender)
phenotype$Prednison <- as.factor(phenotype$Prednison)
phenotype$SmokingBirthMom <- as.factor(phenotype$SmokingBirthMom)
phenotype$SmokingPregnacy <- as.factor(phenotype$SmokingPregnacy)
phenotype$SmokingInsideBirth <- as.factor(phenotype$SmokingInsideBirth)
for (i in 1:n) {
  x <- targets[i,]
  y <- phenotype[phenotype$Sample==targets[i,1],]
  phenotype1[i,] <- c(x[8], x[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11], y[12])
}
# plates are x[3] slides are x[6]
M_matrix <- M.val[,!(duplicated(targets$Sample_Name))]
phenotype <- phenotype1[phenotype1$Basename %in% colnames(M_matrix),]
phenotype$Sample <- NULL
rownames(phenotype) <- phenotype$Basename
phenotype$Basename <- NULL
remove(M.val)
phenotype$Wheeze[phenotype$Wheeze==2] <- 0
phenotype$Asthma[phenotype$Asthma==2] <- 0
phenotype$Prevention[phenotype$Prevention==2] <- 0
phenotype$Gender[phenotype$Gender==2] <- 0
phenotype$Prednison[phenotype$Prednison==2] <- 0
phenotype$SmokingBirthMom[phenotype$SmokingBirthMom==2] <- 0
phenotype$SmokingPregnacy[phenotype$SmokingPregnacy==2] <- 0
phenotype$SmokingInsideBirth[phenotype$SmokingInsideBirth==2] <- 0
# construct an adjacency matrix
k <- which(is.na(M_matrix), arr.ind=TRUE)
M_matrix[k] <- rowMedians(M_matrix, na.rm=TRUE)[k[,1]]
M_matrix<- t(M_matrix)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(M_matrix, dataIsExpr=TRUE, powerVector=powers, verbose=5, networkType="signed")

png(file=paste(loc.wcna, "SFT_wgcna_u.png", sep="/"), width=1980, height=1080)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
setwd(loc.wcna)
beta <- 5
net <- blockwiseModules(M_matrix, power=beta, TOMType="signed", minModuleSize=100, reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE, pamRespectsDendro=FALSE, saveTOMs=TRUE, saveTOMFileBase="maki_tom", verbose=3)
setwd(loc)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
png(file=paste(loc.wcna, "dendro_col_wcna_u.png", sep="/"), width=1980, height=1080)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
methTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, methTree,
file = paste(loc.wcna, "Maki-networkConstruction-auto.RData", sep="/"))

# Define numbers of genes and samples
nGenes = ncol(M_matrix)
nSamples = nrow(M_matrix)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(M_matrix, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, phenotype, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

png(file=paste(loc.wcna, "Maki_wgcna_trait_u.png", sep="/"), width=1980, height=1080)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(phenotype), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
