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
library(sandwich)
library(lmtest)
#sink("/data/f114798/Logs/log_sva.log")
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# set the locations of important files
loc <- "/data/f114798/"
loc.rgdata <- "Data/RG_data"
loc.comp <- "Complementary_data"
loc.mdata <- "Data/M_Values"
loc.model <- "Data/Models/model_9_LungFunc"
setwd(loc)
load(file=paste(loc.comp, "png_settings.Rdata", sep="/"))
load(file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))
M_matrix<- M.val
phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
phenotype <- data.frame(Sample=phenotype[,2], Gender=phenotype[,6], Wheeze=phenotype[,20], FEV05=phenotype[,26], Age=phenotype[,28], SmokingPregnancy=phenotype[,10], stringsAsFactors=False)
n <- length(targets$Sample_Name)
phenotype1 <- data.frame(Basename=rep(NA, n), Sample=rep(NA, n), Gender=rep("", n), Wheeze=rep(NA, n), FEV05=rep(NA, n), Age=rep(NA, n), Batch=rep(NA, n), SmokingPregnancy=rep(NA, n), stringsAsFactors=FALSE)
for (i in 1:n) {
  x <- targets[i,]
  y <- phenotype[phenotype$Sample==targets[i,1],]
  phenotype1[i,] <- c(x[8], x[1], y[2], y[3], y[4], y[5], x[3], y[6])
}
phenotype1$SmokingPregnancy <- as.factor(phenotype1$SmokingPregnancy)
phenotype1$SmokingPregnancy[phenotype1$SmokingPregnancy==2] <- 0
# plates are x[3] slides are x[6]
phenotype <- phenotype1
double.samples <- phenotype$Basename[which(duplicated(phenotype$Sample))]
phenotype <- phenotype[!duplicated(phenotype$Sample),]

# make a model specific dataframe
PHENO <- phenotype[,c("FEV05", "Batch", "Age", "Gender", "SmokingPregnancy")]
# rename the rows and drop the levels from the dataframe
PHENO <- droplevels(PHENO)
rownames(PHENO)<- phenotype[,"Sample"]
PHENO <- PHENO[!is.na(PHENO$FEV05),]
PHENO <- PHENO[!is.na(PHENO$Age),]
PHENO <- PHENO[!is.na(PHENO$SmokingPregnancy),]
# reinstitute the levels and factors
PHENO$Batch <- as.factor(as.numeric(substring(PHENO$Batch, 7)))
PHENO$Age <- as.numeric(PHENO$Age)
PHENO[PHENO[,1]==2,1] <- 0
PHENO$Gender <- as.factor(PHENO$Gender)

# replace NA in M matrix with the median of the row
k = which(is.na(M_matrix), arr.ind=TRUE)  ## replace missing value with median
M_matrix2<- M_matrix;
M_matrix2[k] = rowMedians(M_matrix2, na.rm=TRUE)[k[,1]]
M_matrix <- M_matrix[,colnames(M_matrix) %in% phenotype$Basename]
colnames(M_matrix) <- phenotype$Sample[phenotype$Basename==colnames(M_matrix)]
M_matrix2 <- M_matrix2[,colnames(M_matrix2) %in% phenotype$Basename]
colnames(M_matrix2) <- phenotype$Sample[phenotype$Basename==colnames(M_matrix2)]
M_matrix <- M_matrix[,colnames(M_matrix) %in% rownames(PHENO)]
M_matrix2 <- M_matrix2[,colnames(M_matrix2) %in% rownames(PHENO)]
# create a model and use sva to compare it against the null model
# prepare two models, one with phenotype and all covariates and the other is without phenotype
mod1<-model.matrix(~.,PHENO)
mod0<-model.matrix(~.,PHENO[,-1])
n.sv <- num.sv(dat=M_matrix2, mod=mod1, method = 'leek') # using method leek to estimate number of SVs
# by setting n.sv=NULL in the function sva, the sva will estimate the number automaticly 
svobj1= sva(M_matrix2,mod1,mod0,n.sv=n.sv)
rm(M_matrix2)
SVs = svobj1$sv
colnames(SVs) <-paste0("sv",1:ncol(SVs))

RLMtest= function(methcol, meth_matrix, Y, X1, X2, X3, X4, SV) {
  
  mod= rlm(as.numeric(Y)~meth_matrix[, methcol]+X1+X2+X3+X4+SV,maxit=200)
  cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
  if (class(cf)=="try-error") {
    bad <- as.numeric(rep(NA, 3))
    names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
    bad
  }
  else{
    cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  }
}
M_matrix<- t(M_matrix)
write.table(data.frame(Pheno=dim(PHENO), Mval=dim(M_matrix)), file=paste(loc.model, "info.txt", sep="/"))
system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), RLMtest, meth_matrix=M_matrix, Y=PHENO$FEV05, X1=PHENO$Batch, X2=PHENO$Age, X3=PHENO$Gender, X4=PHENO$SmokingPregnancy, SV<- SVs, mc.cores=12))

all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")

## M_matrix<-t(M_matrix)
all.results<-all.results[match(colnames(M_matrix),all.results$probeID),]
all.results$N<-colSums(!is.na(M_matrix))

# chr
filename <- "MAKI_trimmed_M_glm_model_9"
p.all<- all.results[,c(2:4)]
rownames(p.all)<- all.results[,1]
colnames(p.all)<- c("coef","se","pvalue")
save(p.all,file=paste(loc.model, paste0(filename,".Rdata"), sep="/"))

rownames(all.results)<- as.character(all.results$probeID)
library(GenABEL)
library(GWASTools)

pvalue <- all.results[,4]
names(pvalue)<- as.character(all.results[,1])
padjust <- p.adjust(pvalue, "fdr")
png(file=paste(loc.model, paste0("qqplot_",filename,".png"), sep="/"), width=png.w, height=png.h)
qqPlot(pvalue, main="QQ of methylation model FEV05 = Meth + Batch + Covariables\n + Proxy Variables", cex.lab=font.mp, cex.axis=font.mp, cex.main=font.mp, cex.sub=font.mp)
t<-estlambda(pvalue, method="median",plot=F)
t<- t[[1]]
lamda<- round(t,digit=3)
text (2,5, paste0("lambda=",lamda), cex=font.mp)
dev.off()

num.bonfer<- length(which(p.adjust(pvalue,method="bonferroni")<0.05))
num.cpg <- max(100,num.bonfer)
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
#qqman
Ci<- apply(all.results, 1, function(x)  which(annot$Name==x[1]))
BP <- as.numeric(annot$pos[Ci])
P <- as.numeric(all.results[,4])
CHR <- as.numeric(substring(annot$chr[Ci], 4))
cpg.all <- data.frame(BP, CHR, P)
colnames(cpg.all) <- c("BP", "CHR", "P")
rownames(cpg.all) <- rownames(all.results)
save(cpg.all, file=paste(loc.model, paste0("cpg_all_", filename, ".Rdata"), sep="/"))
png(file=paste(loc.model, paste0("manhattan_plot_", filename, ".png"), sep="/"))
manhattan(cpg.all, col=c("orange1", "azure3"))
dev.off()
