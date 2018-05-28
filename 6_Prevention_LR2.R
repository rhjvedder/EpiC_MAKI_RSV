#!/usr/bin/env Rscript
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(sva)
library(MASS)
library(compare)
library(tableone)
library(matrixStats)
library(plyr)
library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#sink("/data/f114798/Logs/log_sva.log")
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# set the locations of important files
loc <- "/data/f114798/"
loc.rgdata <- "Data/RG_data"
loc.comp <- "Complementary_data"
loc.mdata <- "Data/M_Values"
loc.model <- "Data/Models/model_4_Prevention"
setwd(loc)
load(file=paste(loc.comp, "png_settings.Rdata", sep="/"))
load(file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))
M_matrix<- M.val
phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
phenotype <- data.frame(Sample=phenotype[,2], Gender=phenotype[,6], Wheeze=phenotype[,20], Prevention=phenotype[,27], Age=phenotype[,28], stringsAsFactors=False)
n <- length(targets$Sample_Name)
phenotype1 <- data.frame(Basename=rep(NA, n), Sample=rep(NA, n), Gender=rep("", n), Wheeze=rep(NA, n), Prevention=rep(NA, n), Age=rep(NA, n), Batch=rep(NA, n), stringsAsFactors=FALSE)
for (i in 1:n) {
  x <- targets[i,]
  y <- phenotype[phenotype$Sample==targets[i,1],]
  phenotype1[i,] <- c(x[8], x[1], y[2], y[3], if(y[4]=="Palivizumab") {1} else {0}, y[5], x[3])
}
# plates are x[3] slides are x[6]
phenotype <- phenotype1
M_matrix <- M_matrix
double.samples <- phenotype$Basename[which(duplicated(phenotype$Sample))]
phenotype <- phenotype[!duplicated(phenotype$Sample),]

# make a model specific dataframe
PHENO <- phenotype[,c("Prevention", "Batch")]
# rename the rows and drop the levels from the dataframe
PHENO <- droplevels(PHENO) ; rownames(PHENO)<- phenotype[,"Sample"]
PHENO <- PHENO[!is.na(PHENO$Prevention),]
# reinstitute the levels and factors
batch <- as.factor(PHENO$Batch) ## change coord form (1,2,3) to (1,2)
PHENO$Batch <- as.factor(PHENO$Batch)
PHENO[PHENO[,1]==2,1] <- 0

# replace NA in M matrix with the median of the row
k = which(is.na(M_matrix), arr.ind=TRUE)
M_matrix[k] = rowMedians(M_matrix, na.rm=TRUE)[k[,1]]
M_matrix <- M_matrix[,colnames(M_matrix) %in% phenotype$Basename]
colnames(M_matrix) <- phenotype$Sample[phenotype$Basename==colnames(M_matrix)]

# create a model and use sva to compare it against the null model

GLMtest <- function(methcol, meth_matrix,Y, X1) {
  mod = glm(Y~meth_matrix[, methcol] +X1,family = "binomial")
  cf = summary(mod)$coefficients
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]  
}
M_matrix<- t(M_matrix)
M_matrix <- M_matrix[rownames(M_matrix) %in% rownames(PHENO),]
system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), GLMtest, meth_matrix=M_matrix, Y=PHENO$Prevention, X1=PHENO$Batch, mc.cores=12))

all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")

## M_matrix<-t(M_matrix)
all.results<-all.results[match(colnames(M_matrix),all.results$probeID),]
all.results$N<-colSums(!is.na(M_matrix))

# chr
filename <- "MAKI_trimmed_M_glm_model_4"
p.all<- all.results[,c(2:4)]
rownames(p.all)<- all.results[,1]
colnames(p.all)<- c("coef","se","pvalue")
save(p.all,file=paste(loc.model, paste0(filename,".Rdata"), sep="/"))

rownames(all.results)<- as.character(all.results$probeID)
library(GenABEL)
library(GWASTools)

pvalue<- all.results[,4]
names(pvalue)<- as.character(all.results[,1])
padjust<- p.adjust(pvalue,"fdr")
png(file=paste(loc.model, paste0("qqplot_",filename,".png"), sep="/"), width=png.w, height=png.h)
qqPlot(pvalue, main="QQ of methylation model Intervention = Meth + Batch", cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp)
t<-estlambda(pvalue, method="median",plot=F)
t<- t[[1]]
lamda<- round(t,digit=3)
text (2,5, paste0("lambda=",lamda))
dev.off()

num.bonfer<- length(which(p.adjust(pvalue,method="bonferroni")<0.05))
num.cpg <- max(20,num.bonfer)
cpgnamei<- names(sort(pvalue))[1:num.cpg]
Gindex<- apply(as.matrix(cpgnamei), 1, function(x)  which(annot$Name==x))
genelist<- as.character(annot[Gindex,22])
chr <- as.character(annot$chr[Gindex])
group<- as.character(annot$UCSC_RefGene_Group[Gindex])
cpg.p<- as.character(annot$Relation_to_Island[Gindex])
cpg <- cbind(cpgnamei,all.results[cpgnamei,2], all.results[cpgnamei,3],pvalue[cpgnamei],padjust[cpgnamei], genelist,chr,group,cpg.p)
colnames(cpg)<- c( "Cpgsite", "Beta", "SE", "Pvalue", "P_adjust", "Genelist", "Chr", "Group", "Cpg")
write.csv(cpg, file=paste(loc.model, paste0("cpg_", filename, ".csv"), sep="/"))
save(cpg, file=paste(loc.model, paste0("cpg_", filename, ".Rdata"), sep="/"))


