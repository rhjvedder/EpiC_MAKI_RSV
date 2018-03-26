#!/usr/bin/env Rscript
library(RColorBrewer)
library(pheatmap)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/data/f114798/Logs/log_norm.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/data/f114798/Data"
loc.out.norm <- paste(loc.out, "Norm", sep="/")
norma = "Funnorm"

# loading data
load(file = paste(loc.out, "RG_Channel_Set.Rdata", sep = "/"))
load(file = paste(loc.out, "epic_maki_filtered_gm_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
Phenotype <- read.csv("/data/f114798/Rscripts/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")
targets$Sample_Plate <- trimws(targets$Sample_Plate, "r")
# do a selection on the original rg set with the samples and probes from the gm set
# requires a selection on the original probe names
# rg.set.flt <- rg.set[rownames(rg.set) %in% rownames(gm.set),colnames(rg.set) %in% colnames(gm.set)]
rg.set.flt <- rg.set

norms <- c("Funnorm", "ssNoob", "Quantile")
norms.sets <- sapply(norms, function(norma) {
  # normalize based on the supplied normalizing type
  if(norma == "Funnorm") {
    m.set.sq <- preprocessFunnorm(rg.set.flt)
  } else if(norma == "ssNoob") {
    m.set.sq <- preprocessNoob(rg.set.flt, dyeMethod = "single")
    m.set.sq <- mapToGenome(m.set.sq)
  } else {
    m.set.sq <- preprocessQuantile(rg.set.flt)
  }
  
  save(m.set.sq, file = paste(loc.out, paste("epic_maki_", norma, "_meth_set.Rdata"), sep = "/"))
  # plots showing normalizing before and after, consider making seperate plots for the sample plates
  for (plate in levels(factor(targets$Sample_Plate))) {
    png(paste(loc.out.norm, paste(norma, "_", plate, ".png"), sep = "/"))
    par(mfrow=c(1,2))
    densityPlot(getBeta(rg.set[,rg.set$Basename %in% paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")]),main="Raw", legend=FALSE)
    legend("top", legend = plate,
           text.col=brewer.pal(8,"Dark2"))
    densityPlot(getBeta(m.set.sq[,m.set.sq$Basename %in% paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")]), main="Normalized", legend=FALSE)
    legend("top", legend = plate,
           text.col=brewer.pal(8,"Dark2"))
    print("---------------------------------------------------------------------")
    dev.off()
  }
  return(m.set.sq)
})

print("---------------------------------------------------------------------")

# make a heatmap of the before and after normalizing

nbr <- 1
for (norm.set in norms.sets) {
  norms.betas <- getBeta(norm.set[,norm.set$Sample_Name=="Control_M"])
  plates.order <- sapply(levels(factor(targets$Sample_Plate)), function(plate) {
    basenames <- colnames(norms.betas)[colnames(norms.betas) %in% paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")]
    basenames[order(basenames)]
    basenames
  })
  plates.order <- c(plates.order[1], plates.order[2], plates.order[3], plates.order[4])
  cormat <- cor(norms.betas[,match(plates.order, colnames(norms.betas))])
  png(file=paste(loc.out.norm, paste(norms[nbr], "_heatmap_cor.png"), sep="/"))
  pheatmap(cormat, cluster_rows=FALSE, cluster_cols=FALSE)
  dev.off()
  nbr <- nbr + 1
}

nbr <- 1
png(file=paste(loc.out.norm, "boxplot_all.png", sep="/"))
par(mfrow=c(1,3))
for (norm.set in norms.sets) {
    norms.betas <- getBeta(norm.set[,norm.set$Sample_Name=="Control_M"])
  plates.order <- sapply(levels(factor(targets$Sample_Plate)), function(plate) {
    basenames <- colnames(norms.betas)[colnames(norms.betas) %in% paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")]
    basenames[order(basenames)]
    basenames
  })
  plates.order <- c(plates.order[1], plates.order[2], plates.order[3], plates.order[4])
  cormat <- cor(norms.betas[,match(plates.order, colnames(norms.betas))])
  boxplot(cormat, main=norms[nbr])
  nbr <- nbr + 1
}
dev.off()

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
