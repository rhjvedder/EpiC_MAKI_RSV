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
loc.model <- "Data/Models/model_3_Asthma"
setwd(loc)

load(file=paste(loc.mdata, "MAKI_trimmed_M.Rdata", sep="/"))
M_matrix<- M.val
phenotype <- read.csv(paste(loc.comp, "Phenotype_data_Maki.csv", sep = "/"))
phenotype <- data.frame(Sample=phenotype[,2], Gender=phenotype[,6], Wheeze=phenotype[,20], Asthma=phenotype[,21], Age=phenotype[,28], stringsAsFactors=False)
n <- length(targets$Sample_Name)
phenotype1 <- data.frame(Basename=rep(NA, n), Sample=rep(NA, n), Gender=rep("", n), Wheeze=rep(NA, n), Asthma=rep(NA, n), Age=rep(NA, n), Batch=rep(NA, n), stringsAsFactors=FALSE)
for (i in 1:n) {
  x <- targets[i,]
  y <- phenotype[phenotype$Sample==targets[i,1],]
  phenotype1[i,] <- c(x[8], x[1], y[2], y[3], y[4], y[5], x[3])
}
# plates are x[3] slides are x[6]
phenotype <- phenotype1
M_matrix <- M_matrix
double.samples <- phenotype$Basename[which(duplicated(phenotype$Sample))]
phenotype <- phenotype[!duplicated(phenotype$Sample),]

# make a model specific dataframe
PHENO <- phenotype[,c("Asthma", "Batch", "Age", "Gender")]
# rename the rows and drop the levels from the dataframe
PHENO <- droplevels(PHENO) ; rownames(PHENO)<- phenotype[,"Sample"]
PHENO <- PHENO[!is.na(PHENO$Asthma),]
PHENO <- PHENO[!is.na(PHENO$Age),]
# reinstitute the levels and factors
PHENO$Batch <- as.factor(as.numeric(substring(PHENO$Batch, 7)))
PHENO$Age <- round(PHENO$Age,2)
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
# NA in phenotype data
n.sv <- num.sv(dat=M_matrix, mod=mod1, method = 'leek') # using method leek to estimate number of SVs
# by setting n.sv=NULL in the function sva, the sva will estimate the number automaticly 
svobj1= sva(M_matrix2, mod1, mod0, n.sv=n.sv)
rm(M_matrix2)
SVs = svobj1$sv
colnames(SVs) <-paste0("sv",1:ncol(SVs))

GLMtest <- function(methcol, meth_matrix,Y, X1, X2, X3, SV) {
  mod = glm(Y~meth_matrix[, methcol] +X1+X2+X3+SV,family = "binomial")
  cf = summary(mod)$coefficients
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]  
}
M_matrix<- t(M_matrix)
system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), GLMtest, meth_matrix=M_matrix, Y=PHENO$Asthma, X1=PHENO$Batch, X2=PHENO$Age, X3=PHENO$Gender, SV<- SVs, mc.cores=12))

all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")

## M_matrix<-t(M_matrix)
all.results<-all.results[match(colnames(M_matrix),all.results$probeID),]
all.results$N<-colSums(!is.na(M_matrix))

# chr
filename <- "MAKI_trimmed_M_glm_model_3"
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
png(file=paste(loc.model, paste0("qqplot_",filename,".png"), sep="/"), width=1000, height=1000)
qqPlot(pvalue, main="QQ of methylation model Asthma = Meth + Batch + Covariables\n + Proxy Variables", cex=2)
t <- estlambda(pvalue, method="median",plot=F)
t <- t[[1]]
lamda <- round(t,digit=3)
text(4,1, paste0("lambda=", lamda))
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
