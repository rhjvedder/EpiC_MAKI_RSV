#!/usr/bin/env Rscript
library(minfi)
library(RColorBrewer)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

sink("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Logs/log_2.log")
start.time <- Sys.time()

# example data EpiC 850K
loc.out <- "/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Data"
norma = "Funnorm"

# loading data
load(file = paste(loc.out, "epic_maki_rg_set.Rdata", sep = "/"))
load(file = paste(loc.out, "epic_maki_filtered_rg_set.Rdata", sep = "/"))
manifest <- getManifest(rg.set)
print("manifest")
manifest
Phenotype <- read.csv("/groups/umcg-griac/tmp03/projects/umcg-rhjvedder/Rscript/Phenotype_data_Maki.csv")
print("---------------------------------------------------------------------")

# normalize based on the supplied normalizing type
if(norma == "Funnorm") {
  m.set.sq <- preprocessFunnorm(rg.set.flt)
} else if(norma == "ssNoob") {
  m.set.sq <- preprocessNoob(rg.set.flt, dyeMethod = "single")
  m.set.sq <- mapToGenome(m.set.sq)
} else {
  m.set.sq <- preprocessQuantile(rg.set.flt)
}

save(m.set.sq, file = paste(loc.out, "epic_maki_filtered_meth_set.Rdata", sep = "/"))

# plots showing normalizing before and after
png(paste(loc.out, "normalizing.png", sep = "/"))
par(mfrow=c(1,2))
densityPlot(rg.set, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(m.set.sq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
print("---------------------------------------------------------------------")
dev.off()
print("---------------------------------------------------------------------")

# time the run for later runnin' purposses
print("script took:")
Sys.time() - start.time
# stop logging
sink()
