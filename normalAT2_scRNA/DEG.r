#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell=args[1]
case=args[2]
ctrl=args[3]
wd=args[4]

print('this is a script of MAST filter fix model save to')
print(wd)
print('to identify DEG between')
print(paste(case, 'vs',ctrl, 'in', cell))

library(Seurat)
library(SingleCellExperiment)
library(parallel)
options(mc.cores = 16)
library(stringr)

dat = readRDS('lung_scRNA_normal_AT2_aged_young.rds')

sce <-SingleCellExperiment(assays=list(counts=seurat[['RNA']]@counts), colData=seurat@meta.data)

#normlization TPM
library(scater)
TPM=calculateTPM(sce, exprs_values="counts")
logTPM<-log2(TPM+1)

fData <- data.frame(names=rownames(logTPM))
rownames(fData) <- rownames(logTPM)
cData <- seurat[[]] #data.frame(both=seurat$both, age=seurat$age, condition=seurat$condition, ID=seurat$ID, cluster=seurat$seurat_clusters, sex=seurat$sex, batch=seurat$batch)
rownames(cData) <- colnames(logTPM)

library(MAST)
scaRaw<- FromMatrix(as.matrix(logTPM), cData, fData)

freq_expressed <- 0.1
#FCTHRESHOLD <- log2(1.5)

expressed_genes <- freq(scaRaw) > freq_expressed #genes that are found in at least 0.2 of the sample (since we won⚌~@~Yt have any power to conclude much about less frequent transcripts).
scaRaw <- scaRaw[expressed_genes,]
scaRaw

colData(scaRaw)$cngeneson<-scale(colSums(assay(scaRaw)>0))
colData(scaRaw)$both=factor(colData(scaRaw)$both, levels=c(ctrl,case))

library(reshape2)
library(data.table)
library(NMF)


#testing only the condition, not the replicates
#to interprate the coefficients more easily, we set reference value to Y-SPCCre
#celltype <- relevel(celltype,"AT2") #set reference to AT2
#zlmCond<-zlm(~celltype+cngeneson+sex+both+(1|ID)+(1|batch), scaRaw, method='glmer', ebayes=FALSE) #(1|celltype:ID) interaction
#summary(zlm(~ celltype + (1|Batch), sca = scaRaw, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable %>%
#filter(primerid=="celltypeAT2S")
#zlmCond<-zlm(~celltype+cngeneson+(1|ID), scaRaw, method='glmer', ebayes=FALSE)
zlmCond<-zlm(~both+cngeneson+batch, scaRaw)
setwd(wd)
file=paste0('MAST_filter_fix_', cell,'_',case, ctrl, "_zlmCond.RData")
#save(zlmCond, file=file) #filtre gene name is DEG_MAST_Ased_vs_Ysed_AT2_zlmCond.RData

#zlmCond<-zlm(~age+sex+(1|ID)+(1|batch), scaRaw, method='glmer', ebayes=FALSE, strictConvergence = FALSE)

contrast=paste0('both',case)
summaryCond <- summary(zlmCond, doLRT=paste0('both',case))

file=paste0('MAST_filter_fix_', cell,'_',case, ctrl, '_summaryCond.RData')
#save(summaryCond, file=file)

summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast==paste0('both',case) & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==paste0('both',case) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(scaRaw)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig<-as.data.frame(fcHurdleSig)
#fcHurdleSig<-fcHurdleSig[order(fcHurdleSig$logFC, decreasing=T),]

library(xlsx)
file=paste0('MAST_filter_fix_', cell,'_',case, ctrl, '.xlsx')
write.xlsx(fcHurdleSig, file=file)
