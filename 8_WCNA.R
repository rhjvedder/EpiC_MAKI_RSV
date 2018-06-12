#!/usr/bin/env Rscript
library(WGCNA)
library(plyr)
library(parallel)
library(matrixStats)
options(stringsAsFactors = FALSE)
# set the locations of important files
allowWGCNAThreads(nThreads = 12)
loc <- "/data/f114798/"
loc.comp <- "Complementary_data"
loc.mdata <- "Data/M_Values"
loc.wcna <- "Data/Wcna/subset_var.4"
setwd(loc)
load(file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))
load(file=paste(loc.comp, "png_settings.Rdata", sep="/"))
M_matrix <- M.val
rm(M.val)
phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
phenotype <- data.frame(Sample=phenotype[,2], Wheeze=phenotype[,20], Asthma=phenotype[,21], Intervention=phenotype[,27], FEV05=phenotype[,26], Gender=phenotype[,6], Age=phenotype[,28], IgE=phenotype[,19], Prednison=phenotype[,22],  SmokingPregnacy=phenotype[,10], stringsAsFactors=FALSE)
n <- length(targets$Sample_Name)
phenotype1 <- data.frame(Basename=rep(NA, n), Sample=rep(NA, n), Wheeze=rep(NA, n), Asthma=rep(NA, n), Intervention=rep(NA, n), FEV05=rep(NA, n), Gender=rep(NA, n), Age=rep(NA, n), IgE=rep(NA, n), Prednison=rep(NA, n),  SmokingPregnacy=rep(NA, n), stringsAsFactors=FALSE)

phenotype$Wheeze <- as.factor(phenotype$Wheeze)
phenotype$Asthma <- as.factor(phenotype$Asthma)
phenotype$Intervention <- as.factor(phenotype$Intervention)
phenotype$Gender <- as.factor(phenotype$Gender)
phenotype$Prednison <- as.factor(phenotype$Prednison)
phenotype$SmokingPregnacy <- as.factor(phenotype$SmokingPregnacy)
for (i in 1:n) {
  x <- targets[i,]
  y <- phenotype[phenotype$Sample==targets[i,1],]
  phenotype1[i,] <- c(x[8], x[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10])
}
# plates are x[3] slides are x[6]
phenotype <- phenotype1[phenotype1$Basename %in% colnames(M_matrix),]
phenotype$Sample <- NULL
rownames(phenotype) <- phenotype$Basename
phenotype$Basename <- NULL
phenotype$Wheeze[phenotype$Wheeze==2] <- 0
phenotype$Asthma[phenotype$Asthma==2] <- 0
phenotype$Intervention[phenotype$Intervention==2] <- 0
phenotype$Gender[phenotype$Gender==2] <- 0
phenotype$Prednison[phenotype$Prednison==2] <- 0
phenotype$SmokingPregnacy[phenotype$SmokingPregnacy==2] <- 0
# construct an adjacency matrix
k <- which(is.na(M_matrix), arr.ind=TRUE)
M_matrix[k] <- rowMedians(M_matrix, na.rm=TRUE)[k[,1]]
M_matrix<- t(M_matrix)
# select only probes with a variance higher than 1
var_n = 0.35
M_matrix <- as.matrix(M_matrix[,which(apply(M_matrix,2,var)> var_n )])
probes.n <- nrow(M_matrix)

# Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# sft = pickSoftThreshold(M_matrix, dataIsExpr=TRUE, powerVector=powers, verbose=5, networkType="signed")

# png(file=paste(loc.wcna, "SFT_wgcna_u.png", sep="/"), width=1980, height=1080)
# par(mfrow = c(1,2))
# font.mp = 2
# cex1 = 1.5
# Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"), cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp);
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"), cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp)
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()

setwd(loc.wcna)
beta <- 12
net <- blockwiseModules(M_matrix, power=beta, TOMType="signed", minModuleSize=10, reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE, pamRespectsDendro=FALSE, saveTOMs=FALSE, verbose=3, maxBlockSize=60000)
setwd(loc)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
png(file=paste(loc.wcna, "dendro_col_wcna_u.png", sep="/"), width=1980, height=1080)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp)
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
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(phenotype), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"), cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp)
dev.off()

print(var_n)
print(beta)
print(probes.n)
