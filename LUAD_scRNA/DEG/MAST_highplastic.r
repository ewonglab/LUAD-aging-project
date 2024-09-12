#use similar model as simon to put cell type, sex,batch and age ie age + cl_1st + sex + batch
#except filter and add cellular detection rate
library(Seurat)
library(SingleCellExperiment)
setwd('/g/data/zk16/qing/tuomas/simon/')

all=readRDS('202111251525_X_LUAD_aging_pp_sce.rds')

meta=read.csv('202111251525_X_LUAD_aging_KP_tumor_subsampleUAD_aging_KP_tumor_subsample_meta.csv', row.names=1)
#meta=read.csv('202111251525_X_LUAD_aging_pp_meta.csv')#, row.names=1)

library(dplyr)

#meta1=left_join(meta1, meta, by='X')

#length(intersect(meta1$X, meta$X))  #complete overlap as expected

d=all[,rownames(meta)]
rm(all)

d$cl_epithelial_tumor_pca=meta$cl_epithelial_tumor_pca #add cell type
sub_meta=subset(meta,  cl_epithelial_tumor_pca=='High plasticity cell state' & tumor_stage != 'healthy')
 sub=d[,match(rownames(sub_meta), colnames(d))]
#  count_raw = assay(sub, "raw")#as.matrix(sub@assays$RNA@counts)
#count_raw=assay(d, "raw")
#meta=colData(d)
#subset the cell type and condition for comparisons  
 #seurat=subset(dat, celltype ==cell &  age %in% c(case,ctrl))

# create a single cell experiment object for downstream processing
#sce <-SingleCellExperiment(assays=list(counts=count_raw, colData=meta))

#normlization TPM
library(scater)
TPM=calculateTPM(sub, exprs_values="raw")
logTPM<-log2(TPM+1)

fData <- data.frame(names=rownames(logTPM))
rownames(fData) <- rownames(logTPM)
cData <- colData(sub)#seurat[[]]
rownames(cData) <- colnames(logTPM)

#filter first
#rowSums(logTPM>0)/

library(MAST)
scaRaw<- FromMatrix(as.matrix(logTPM), cData, fData)

freq_expressed <- 0.1
#FCTHRESHOLD <- log2(1.5)

expressed_genes <- freq(scaRaw) > freq_expressed #genes that are found in at least 0.1 of the sample (since we won⚌~@~Yt have any power to conclude much about less frequent transcripts).
scaRaw <- scaRaw[expressed_genes,]
scaRaw

case='aged'
ctrl='young'
colData(scaRaw)$cngeneson<-scale(colSums(assay(scaRaw)>0)) #cellular detection rate
colData(scaRaw)$batch=factor(colData(scaRaw)$batch)
colData(scaRaw)$sex=factor(colData(scaRaw)$sex)
colData(scaRaw)$cl_epithelial_tumor_pca=factor(colData(scaRaw)$cl_epithelial_tumor_pca)
colData(scaRaw)$age=factor(colData(scaRaw)$age, levels=c(ctrl,case))


######### THIS PART MUST BE DONE WITH MULTIPLE CPUs.
library(reshape2)
library(data.table)
library(NMF)

zlmCond<-zlm(~age+cngeneson+batch+sex, scaRaw) #only one cell type so no need to include cl_epithelial_tumor_pca

contrast=paste0('age',case)
summaryCond <- summary(zlmCond, doLRT=paste0('age',case))

summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast==paste0('age',case) & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==paste0('age',case) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(scaRaw)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig<-as.data.frame(fcHurdleSig)
#write.csv(fcHurdleSig, file='/g/data/zk16/qing/tuomas/simon/MAST/MAST_highplastic_17wk_agedvsyoung.csv')

library(xlsx)
write.xlsx(fcHurdleSig, file="/g/data/zk16/qing/tuomas/simon/MAST/MAST_highplastic_agedvsyoung.xlsx")


