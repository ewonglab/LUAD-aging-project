
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  print("the batch number has been provided, start to run Rscript {batchtime}")
}

#modified on emily's script 2101
# library
library(DropletUtils)
library(Seurat)
library(ggplot2)
library(scran)

#rundate='XZ-2284_20230405_Hyperoxia_day1'
rundate= args[1]
setwd('/g/data/zk16/projects/tammela/R3_shNupr1/')
fname <- file.path( paste0( rundate,'/outs/raw_feature_bc_matrix'))
sce <- read10xCounts(fname, col.names=TRUE) #6794880
hash=tail(rownames(sce), 10)
sce= sce[!rownames(sce) %in% hash, ]

#--- gene-annotation ---#
library(scater)
rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$ID, rowData(sce)$Symbol)
library(EnsDb.Mmusculus.v79)
location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(sce)$ID, #Unable to map 12284 of 54855 requested IDs.
    column="SEQNAME", keytype="GENEID")

#--- cell-detection test ---#
set.seed(100)
limit <- 100
all.out <- emptyDrops(counts(sce), lower=limit, test.ambient=TRUE)
setwd('/g/data/zk16/qing/tuomas/shNupr1_scRNA/03emptydrop/')
pdf(paste0( rundate,'_pvalue_distribution.pdf'))
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
    xlab="P-value", main="", col="grey80")
dev.off()


#--- cell-detection ---#
set.seed(100)
e.out <- emptyDrops(counts(sce))
sce <- sce[,which(e.out$FDR <= 0.001)]
summary(e.out$FDR <= 0.001)
#   Mode   FALSE    TRUE    NA's
#logical   53158   15591 6726131

is.cell <- e.out$FDR <= 0.001
nonemptycell = colnames(sce) #15591
write.table(nonemptycell, file=paste0(rundate,"_nonemptycell.txt"))


#--The distribution of total counts--#
library(DropletUtils)
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)

pdf(paste0( rundate, '_UMIdistribution.pdf'))
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"),
        col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

dev.off()


#Note - if you only want the cell IDs to remove you can just save these cellIDs
#https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
pdfname=paste0("emptydrops_plot_", rundate, ".pdf")
pdf(pdfname)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

# export to seurat remove scublet and hash filtered cells compare to original dataset
setwd('/g/data/zk16/qing/tuomas/shNupr1_scRNA/01dehash/')
rdsname=paste0(rundate, ".rds")
R1=readRDS(rdsname) # where cells filtered based on hashing information
head(colnames(R1))

length(colnames(R1)) #6771

#remove scrublet cells
setwd('/g/data/zk16/qing/tuomas/shNupr1_scRNA/02scrublet/')
doub = read.delim('all_doublet_cellID.txt',header=F, sep=' ') 
R1 <- R1[,! colnames(R1) %in% doub$V3] 

print("number of cell left after removing scrublet doublet cells is")
length(colnames(R1)) 

#remove empty cells
R1 <- R1[, colnames(R1) %in% nonemptycell]

print("number of cell left after removing empty cells is")
length(colnames(R1)) #6211

#QC
library(tibble)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
R1[["percent.mt"]] <- PercentageFeatureSet(R1, pattern = "^mt-")
summary(R1[["percent.mt"]] )

# Percent of mitochondrial ribosomal https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html
ribo.genes <- grep(pattern = "(^Rpl|^Rps|^Mrp)", x = rownames(R1[['RNA']]), value = TRUE)
percent.ribo <- Matrix::colSums(R1[['RNA']][ribo.genes, ])/Matrix::colSums(R1[['RNA']])
#Raw.data <- AddMetaData(object = Raw.data, metadata = percent.ribo, col.name = "percent.ribo")

R1[["percent.ribo"]] <- percent.ribo
summary(R1@meta.data$percent.ribo)


# Visualize QC metrics as a violin plot
setwd('/g/data/zk16/qing/tuomas/shNupr1_scRNA/04QC/')

pdf(paste0(rundate,"QC_before_filter.pdf"))
VlnPlot(R1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
dev.off()

#before  QC filter
rdsname=paste0(rundate,"_noempty_nodoublet.rds")
saveRDS(R1, file=rdsname)


#filter

R1 <- subset(R1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

pdf(paste0(rundate,"QC_after_filter.pdf"))
VlnPlot(R1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
dev.off()


#-------------scran normalization
#https://bioconductor.org/books/release/OSCA/normalization.html


##----standar seurat pipeline
R1 <- FindVariableFeatures(R1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(R1)
R1 <- ScaleData(R1, features = all.genes)
R1 <- RunPCA(R1, features = VariableFeatures(object = R1))


R1 <- FindNeighbors(R1, dims = 1:10)
R1 <- FindClusters(R1, resolution = 0.5)
R1 <- RunUMAP(R1, dims = 1:10)

features=c('Axin2', 'Tm4sf1', 'Cp', 'Hhip', 'Wnt5a', 'Sftpb', 'Sftpc', 'Lyz1', 'Lyz2', 'H2-Ab1')
pdfname=paste0(rundate, "_noempty_nodoublet_feature_umap.pdf")
pdf(pdfname)
DimPlot(R1, reduction = "umap")
DimPlot(R1, reduction = "umap", group.by='Age')
DimPlot(R1, reduction = "umap", group.by='shRNA.group')
DimPlot(R1, reduction = "umap", group.by='Sex')
dev.off()

pdfname=paste0(rundate, "_noempty_nodoublet_feature_violin.pdf")

pdf(pdfname, width=16)
VlnPlot(R1, features = features, pt.size = 0)
dev.off()

rdsname=paste0(rundate,"_noempty_nodoublet_qc.rds")
saveRDS(R1, file=rdsname)

