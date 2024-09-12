library(openxlsx)
de=read.xlsx('MAST_filter_fix_AT2_aged_sedyoung_sed.xlsx', rowNames = T) #wtAT2 DEG
colnames(de)[1]='gene' 
rownames(de)=de$gene

setwd('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/Simon/LUADversion2/MAST/celltype')
fn=list.files(pattern = '.xlsx') # LUAD DEG per cell type
fn
f=lapply(fn, read.xlsx)
names(f)=sapply(strsplit(fn,"_"), `[`, 2)
f$AT2=de #combine
f=f[c('AT2', sapply(strsplit(fn,"_"), `[`, 2))]
lapply(f, head)

#combine all
library(dplyr)
library(tidyverse)
df=f %>% reduce(inner_join, by = "names")
rownames(df)=df$gene
head(df)

subplt=df %>% subset(fdr.x<0.1) #select the fdr<0.1 in AT2

subplt=subplt[,grep('^coef*',colnames(df))]

#top 25 upregulated , select based on the fdr ranking of AT2 only 
up=subplt[rowSums(subplt>0)==ncol(subplt),] #select the ones with all coef> 0 
up=up[order(abs(rowSums(up)),decreasing=T),]
 
#reorder to make it from small log2fc to large log2fc
up=up[1:25,]
up=up[order(abs(rowSums(up)),decreasing=F),]
rownames(up)

#top 25 downregulated
down=subplt[rowSums(subplt<0)==ncol(subplt),]  #select the ones with all coef> 0 
down=down[order(abs(rowSums(down)),decreasing=T),]

down=down[1:25,]
rownames(down)

#combine and prepare the plot data frame
plt=rbind( down, up)

colnames(plt)=names(f)
plt=t(plt)

plt=plt[c('AT2', 'AT2like', 'AT1like', 'highplastic', 'Endodermlike', 'Ribosomehigh'),]

rownames(plt)=c('AT2', 'AT2-like', 'Hopx+ intermediate', 'High plasticity state', 'Endoderm-like', 'Ribosome high')

metadata=matrix(nrow=ncol(plt),
                ncol=1)

colnames(metadata)='change'
rownames(metadata)=colnames(plt)

metadata[,'change']=ifelse(rownames(metadata) %in% rownames(up), 'aged' , 'young')
  
metadata


#Set annotation
metadata=as.data.frame(metadata)
ann <- data.frame(metadata$change)
colnames(ann) <- 'change'
colours <- list('change' = c('aged' = 'red', 'young' = 'darkblue'))

colAnn <- HeatmapAnnotation(df = ann,
  which = 'col',
  col = colours,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'))

paletteLength=10
col1 = "darkblue"
    col2 = 'white'
    col3 = 'red'
    
myBreaks <- c(seq(min(plt), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(plt)/paletteLength, max(plt), length.out=floor(paletteLength/2)))
myBreaks

col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( length(myBreaks))
col.heatmap

pdf('heatmap_MASTcoef_most_correlated_DEG_colAnn_rowSumsrank_colorcode_gene_top.pdf', width = 20)
Heatmap(plt, cluster_columns=FALSE, cluster_rows = F, show_row_names = T, show_column_names = T,  
        #top_annotation_height=unit(1,'cm'), 
        top_annotation=colAnn,
         column_names_side = c("top"),
        width = ncol(plt)*unit(5, "mm"), 
    height = nrow(plt)*unit(5, "mm"),
    #column_order=keep[order(de[keep,]$coef, decreasing=T)], 
    border=T, #name = 'Odds Ratio', #column_title =main, #right_annotation = ha,
        col =colorRamp2(myBreaks, col.heatmap))
dev.off()
