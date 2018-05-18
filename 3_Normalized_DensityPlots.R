#!/usr/bin/env Rscript
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

loc <- "/data/f114798/"
loc.rgdata <- "Data/RG_data"
loc.norm <- "Data/Norm"
loc.data <- "Data/meth_set"
setwd(loc)

load(file = paste(loc.rgdata, "RG_Channel_Set.Rdata", sep = "/"))
targets$Sample_Plate <- trimws(targets$Sample_Plate, "r")

norms <- c("Funnorm", "ssNoob", "Quantile")

norms.sets <- sapply(norms, function(norma) {
  load(paste(loc.data, paste("epic_maki_", norma, "_meth_set.Rdata"), sep = "/"))
  return(m.set.sq)
})

getSampleSet <- function(set) {
  return(set[rownames(set) %in% sample(rownames(set), 50000),])
}

rg.set <- preprocessRaw(rg.set)

i <- 1
for (norma in norms) {
  m.set.sq <- norms.sets[[i]]
# plots showing normalizing before and after, consider making seperate plots for the sample plates
  for (plate in levels(factor(targets$Sample_Plate))) {
    type.set <- getProbeType(rg.set)
    samples <- paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")
    raw.set <- getSampleSet(getBeta(rg.set[,rg.set$Basename %in% samples]))
    raw.t1 <- getSampleSet(getBeta(rg.set[type.set=="I", rg.set$Basename %in% samples]))
    raw.t2 <- getSampleSet(getBeta(rg.set[type.set=="II", rg.set$Basename %in% samples]))
    norm.set <- getSampleSet(getBeta(m.set.sq[,m.set.sq$Basename %in% samples]))
    norm.t1 <- getSampleSet(getBeta(m.set.sq[type.set=="I", m.set.sq$Basename %in% samples]))
    norm.t2 <- getSampleSet(getBeta(m.set.sq[type.set=="II", m.set.sq$Basename %in% samples]))
    sets.names <- list("Raw", "Raw T1", "Raw T2", "Normalized", "Norm T1", "Norm T2")
    sets.data <- list(raw.set, raw.t1, raw.t2, norm.set, norm.t1, norm.t2)
    png(paste(loc.norm, paste(norma, "_probeType_", plate, ".png"), sep = "/"))
    par(mfrow=c(2,3))
    for (n in 1:6) {
      densityPlot(sets.data[[n]], main=sets.names[[n]], legend=FALSE)
      legend("top", legend = plate, text.col=brewer.pal(8,"Dark2"))
      n <- n + 1
    }
    dev.off()
  }
  i <- i + 1
}

# make a heatmap of the before and after normalizing
for (norma in norms) {
  norm.set <- norms.sets[norma]
  norms.betas <- getBeta(norm.set[,norm.set$Sample_Name=="Control_M"])
  plates.order <- sapply(levels(factor(targets$Sample_Plate)), function(plate) {
    basenames <- colnames(norms.betas)[colnames(norms.betas) %in% paste(targets$Slide[targets$Sample_Plate==plate], targets$Array[targets$Sample_Plate==plate], sep="_")]
    basenames[order(basenames)]
    basenames
  })
  plates.order <- c(plates.order[1], plates.order[2], plates.order[3], plates.order[4])
  cormat <- cor(norms.betas[,match(plates.order, colnames(norms.betas))])
  save(cormat, file=paste(loc.norm, paste(norma, "_heatmap_cor.png"), sep="/"))
  png(file=paste(loc.norm, paste(norms[nbr], "_heatmap_cor.png"), sep="/"))
  pheatmap(cormat, cluster_rows=FALSE, cluster_cols=FALSE)
  dev.off()
}

cormat <- stack(data.frame(sapply(norms.sets, function(norm.set) {
  norms.betas <- getBeta(norm.set[,norm.set$Sample_Name=="Control_M"])
  return(as.vector(cor(norms.betas)))
})))
colnames(cormat) <- norms
png(file=paste(loc.norm, "boxplot_controls.png", sep="/"))
ggplot(dat) + 
  geom_boxplot(data=cormat, aes(x = ind, y = values))
dev.off()
