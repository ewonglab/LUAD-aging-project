library(Seurat)
# Load in the UMI matrix and HTO count matrix


runs= c('IGO_028755', 'IGO_028756', 'IGO_028761')
for (i in 1:length(runs)){
rundate= runs[i]
print(rundate)

setwd('/g/data/zk16/projects/tammela/R2_shDNMT1/')
datadir= paste0( './', rundate, '/outs/filtered_feature_bc_matrix/')

expression_matrix = Read10X(
	data.dir =datadir,
	gene.column = 2,
	unique.features = TRUE,
	strip.suffix = FALSE
	)

	
setwd('/g/data/zk16/qing/tuomas/shDNMT1_scRNA/02scrublet/')
cellname=colnames(expression_matrix$`Gene Expression`)
filename=paste0(rundate,'_doublet_scores.txt')
scrublet_score<-read.table(filename)
head(scrublet_score)
colnames(scrublet_score)="score"
scrublet_score$cell=cellname

library(dplyr)
scrublet_score=scrublet_score %>% filter(score>0.25)
dim(scrublet_score) #273 #243

filename=paste0(rundate,'_doublet_scores_doublet.txt')

write.table(scrublet_score, file=filename, quote=F, row.names=T, col.names=T) 


}


