# installation of bioconductor, minfi and loading of libraries
# source("https://bioconductor.org/biocLite.R")
# install.packages("RColorBrewer")
# biocLite("minfi")
# 1biocLite("minfiDataEPIC")
library(minfi)
library(minfiData)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# example data EpiC 850K


makie_pipe <- function(loc, norma="Quantile") {
  base.dir <- system.file("extdata", package = loc)
  list.files(base.dir)
  targets <- read.metharray.sheet(base.dir)
  rg.set <- read.metharray.exp(targets = targets)
  
  # Annotation
  ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  # QC on signal detection
  det.p <- detectionP(rg.set)
  keep <- colMeans(det.p) < 0.05
  rg.set <- rg.set[,keep]
  targets <- targets[keep,]
  det.p <- det.p[,keep]
  
  # normalize
  if(norma == "Funnorm") {
    m.set.sq <- preprocessFunnorm(rg.set)
  } else if(norma == "ssNoob") {
    m.set.sq <- preprocessNoob(rg.set, dyeMethod = "single")
    m.set.sq <- mapToGenome(m.set.sq)
  } else {
    m.set.sq <- preprocessQuantile(rg.set)
  }

  # visualize before and after normalizing
  par(mfrow=c(1,2))
  densityPlot(rg.set, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
  legend("top", legend = levels(factor(targets$Sample_Group)),
         text.col=brewer.pal(8,"Dark2"))
  densityPlot(getBeta(m.set.sq), sampGroups=targets$Sample_Group,
              main="Normalized", legend=FALSE)
  legend("top", legend = levels(factor(targets$Sample_Group)),
         text.col=brewer.pal(8,"Dark2"))
  
  # ensure probes are in the same order in the m.set.sq and det.p objects
  det.p <- det.p[match(featureNames(m.set.sq),rownames(det.p)),]
  
  # remove any probes that have failed in one or more samples
  keep <- rowSums(det.p < 0.01) == ncol(m.set.sq)
  m.set.sq.flt <- m.set.sq[keep,]
  
  # add sex to samples
  predicted.sex <- getSex(m.set.sq.flt, cutoff = -2)
  plotSex(predicted.sex)
  #m.set.sq.flt <- addSex(m.set.sq.flt, predicted.sex$predictedSex)
  
  # if your data includes males and females, remove probes on the sex chromosomes
  keep <- !(featureNames(m.set.sq.flt) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
  m.set.sq.flt <- m.set.sq.flt[keep,]
  
  # filter cross reactive probes
  reactive.probes <- read.csv(file="1-s2.0-S221359601630071X-mmc1.csv", sep="\t", stringsAsFactors=FALSE)
  keep <- !(featureNames(m.set.sq.flt) %in% reactive.probes$IlmnID)
  m.set.sq.flt <- m.set.sq.flt[keep,]
  
  # filter polymorphic targets
  
  polymorphic.probes <- read.csv(file="1-s2.0-S221359601630071X-mmc2.csv", sep="\t", stringsAsFactors=FALSE)
  keep <- !(featureNames(m.set.sq.flt) %in% polymorphic.probes$IlmnID)
  m.set.sq.flt <- m.set.sq.flt[keep,]
  
  # remove probes with SNPs at CpG site
  m.set.sq.flt <- dropLociWithSnps(m.set.sq.flt, maf=0.05)
  
  # at a filtering step to take out certain SNPs with rsxxxxx
  # load the list probably or search for SNPs with /\
  # remove them
  
  # maybe change the order of filtering and normalizing
  return <- m.set.sq.flt
  
  #je suis f***ed
}

mset.sq.flt <- makie_pipe("minfiDataEPIC", "ssNoob")
