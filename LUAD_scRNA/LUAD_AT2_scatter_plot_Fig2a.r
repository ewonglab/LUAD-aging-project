library(openxlsx)
s=read.xlsx('MAST_tumorcells_nohealthy_agedvsyoung_addribo.xlsx', rowNames = T) 
head(s)
colnames(s)[1]='gene' 

de=read.xlsx('MAST_filter_fix_AT2_aged_sedyoung_sed.xlsx', rowNames = T)
head(de)
colnames(de)[1]='gene' 

library(dplyr)
df=inner_join(s,de,by='gene',suffix=c('.luad', '.normal')) 
head(df)

cor.test(df$coef.exercise, df$coef.luad)

df=as.data.frame(df)
rownames(df)=df$gene

subplt=plt[which(plt$fdr.at2<0.05 & plt$fdr.luad<0.05),]

keep=c('Nupr1', 'Lcn2')
subplt$color=ifelse(subplt$gene %in% keep, "red", "black")
subplt$label=ifelse(subplt$gene %in% keep, subplt$gene, '')

p3=ggscatter(subplt, x='coef.at2', y='coef.luad', label = 'label', 
             color = 'color', palette = c('black', 'red'),
             add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE )+
   xlab('log2FC(aged sedentary/young sedentary) AT2')+
  ylab('log2FC(aged sedentary/young sedentary) LUAD')+
  ggtitle('DEG\nfdr<0.05 in both')+
    stat_cor(method = "pearson", label.x = -2)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
font("legend.title",size = 20, face = "bold")+ #increase size
 font("legend.text", size = 20,face = "bold")+
 font("caption", size = 20, face = "bold")+
 font("xlab", size = 20, face = "bold")+
 font("ylab", size = 20, face = "bold")+
  font("title", size = 20, face = "bold")+
 font("xy.text", size = 20, face = "bold")+
  theme(legend.position = 'none')
p3

pdf('log2FC_luad_tumorcell_nohealthy_ribo_at2_v2_label_sigboth_line_lcn2.pdf')
p3
dev.off()

subplt$color=ifelse(subplt$coef.at2>0 & subplt$coef.luad>0, '++', 
                 ifelse(subplt$coef.at2>0 & subplt$coef.luad<0, '+-',
                        ifelse(subplt$coef.at2<0 & subplt$coef.luad<0, '--', '-+')))
table(subplt$color)

subplt %>% 
  group_by(color) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
