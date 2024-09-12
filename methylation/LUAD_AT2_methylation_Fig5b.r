# read in beta values 
bluad=readRDS('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/methylation/Project_13012_C_LungCancer_Processed/Results/explore/beta_value.rds')

bat2=readRDS('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/methylation/Project_13012_B/Results/explore/beta_value.rds')

#add coordinates
rat2=readRDS('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/methylation/Project_13012_B/Results/explore/beta_value_ranges.rds')

rluad=readRDS('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/methylation/Project_13012_C_LungCancer_Processed/Results/explore/beta_value_ranges.rds')

#read in top 25 up/down regulated DEG name
up=read.table('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/Simon/LUADversion2/MAST/newplot/heatmap_LUAD_normalAT2_DEG_log2FC_most_correlated_DEG_up_in_aged.txt')

down=read.table('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/Simon/LUADversion2/MAST/newplot/heatmap_LUAD_normalAT2_DEG_log2FC_most_correlated_DEG_down_in_aged.txt')

keep=c(rownames(down)[1:25], rownames(up)[25:1]) #select the top 25 genes 

#read in the ensemble gtf files and select the gene coordinates and only keep methylation site that overlap with those regions 
prom=readRDS('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/useful/mm10/mm10.ensGene.gtf_promoters_gr.rds')
prom

#change ensmus gene name to gene symbol
library("org.Mm.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Mm.eg.db, keys = prom$gene_name, keytype = "ENSEMBL", column="SYMBOL")

prom$symbols=symbols
prom

#select interesting genes 
promkeep=prom[which(prom$symbols %in% keep)]

# find sites overlap with promoter
df=findOverlaps(rat2, promkeep)
df

at2=bat2[queryHits(df),] 

at2=as.data.frame(at2)

at2$gene=promkeep$symbols[subjectHits(df)]

at2$young=rowMeans(at2[,c('YM_8_13', 'YM_8_12', 'YF_8_10', 'YM_8_17')])
at2$aged=rowMeans(at2[,c('AF_8_10', 'AM_8_13','AM_8_12','AM_8_17')])

head(at2)

library(dplyr)
plt1=at2 %>% group_by(gene) %>% summarise(meanY=mean(young))
plt2=at2 %>% group_by(gene) %>% summarise(meanA=mean(aged))

plt1

plt=inner_join(plt1, plt2, by='gene')
plt

plt$diff=plt$meanA-plt$meanY
pltat2=plt

df=findOverlaps(rluad, promkeep)
df

luad=bluad[queryHits(df),]
luad=as.data.frame(luad)


luad$gene=promkeep$symbols[subjectHits(df)]

luad$young=rowMeans(luad[,5:9]) #the first 4 are aged and the rest are young 
luad$aged=rowMeans(luad[,1:4])


head(luad)

library(dplyr)
plt1=luad %>% group_by(gene) %>% summarise(meanY=mean(young))
plt2=luad %>% group_by(gene) %>% summarise(meanA=mean(aged))

plt1

plt=inner_join(plt1, plt2, by='gene')
plt

plt$diff=plt$meanA-plt$meanY
pltluad=plt
pltluad

plt=inner_join(pltat2, pltluad, by='gene', suffix=c('.at2','.luad'))

plt$color=ifelse(plt$gene %in% rownames(up), 'aged', 'young')

head(plt)

library(ggpubr)
pdf('C:/Users/q.wang/OneDrive - Victor Chang/Documents/tuomas/methylation/Project_13012_C_LungCancer_Processed/compare/target_beta_value_interesting_genes_comparison.pdf')
ggscatter(plt, x='diff.at2', y='diff.luad',color='color', palette = c('red','blue'),
          add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE )+stat_cor(method = "pearson", label.x = -0.05)+
  xlab('differential methylation value A-Y (AT2)')+
  ylab('differential methylation value A-Y (LUAD)')+
  font("legend.title",size = 20, face = "bold")+ #increase size
 font("legend.text", size = 20,face = "bold")+
 font("caption", size = 20, face = "bold")+
 font("xlab", size = 20, face = "bold")+
 font("ylab", size = 20, face = "bold")+
  font("title", size = 20, face = "bold")+
 font("xy.text", size = 20, face = "bold")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")
dev.off()
