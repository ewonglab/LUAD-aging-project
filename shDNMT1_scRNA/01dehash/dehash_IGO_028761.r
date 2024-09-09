#https://satijalab.org/seurat/archive/v3.0/hashing_vignette.html
library(Seurat)
# Load in the UMI matrix and HTO count matrix
setwd('/g/data/zk16/projects/tammela/R2_shDNMT1/')
expression_matrix = Read10X(
        data.dir ='./IGO_028761/outs/filtered_feature_bc_matrix/',
        gene.column = 2,
        unique.features = TRUE,
        strip.suffix = FALSE
        )


obj = CreateSeuratObject(counts = expression_matrix$`Gene Expression`)
obj[['HTO']] = CreateAssayObject(counts = expression_matrix$`Antibody Capture`)
obj #8294

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")

#filter cells
#log2(as.matrix(obj@assays$HTO@counts)+1) > 5
library(reshape2)
plotdata=melt(log2(as.matrix(obj@assays$HTO@counts)+1))

library(ggplot2)
#ggplot(plotdata, aes(x=value))+geom_histogram()+ facet_wrap(~Var1,scales='free_x')
setwd('/g/data/zk16/qing/tuomas/shDNMT1_scRNA/01dehash')
#ggsave("hto_IGO_028761_filtered_feature_bc_matrix.pdf")

htos=log2(as.matrix(obj@assays$HTO@counts)+1)
thre=c(5.2, 5.5, 5.2,5,
      5, 4, 6,6,
     6, 6 )
for (i in 1:nrow(htos)){
        htos[i,]<-htos[i,]>thre[i]

}

htos<-htos[,colSums(htos) == 1]
dim(htos) #6771

obj=obj[,colnames(htos)]

#get meta info
htos=apply(htos, 2, function(x) ifelse(x==1, T, F))

hashid=rownames(htos)
hash=c()
for (i in 1:ncol(htos)){
        hash[i]=hashid[htos[,i]]
}

hash=as.data.frame(hash)
table(hash$hash)

#read in meta file
library(xlsx)
bc_info=read.xlsx("/g/data/zk16/qing/tuomas/shDNMT1_scRNA/shDNMT1_scRNA_meta.xlsx",sheetIndex=1)
#subset to the batch

library(dplyr)
bc_info=bc_info %>% filter(X10x.submission.ID=="IGO-028761")

bc_info$hash=paste0(bc_info$Hashtag, '-TotalSeqB')

meta=full_join(hash, bc_info, by=c("hash"))
rownames(meta)=colnames(htos) #add cell names before add meta

obj=AddMetaData(
  object = obj,
  metadata = meta,
  col.name = colnames(meta)
)
obj
saveRDS(obj, file="IGO_028761.rds")

