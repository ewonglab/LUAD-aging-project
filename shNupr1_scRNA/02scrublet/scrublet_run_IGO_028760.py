

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import scrublet as scr
import pandas as pd 


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

#Load counts matrix and gene list
input_dir = '/g/data/zk16/projects/tammela/R3_shNupr1/'
rundate= 'IGO_028760'
folder_dir = 'outs/filtered_feature_bc_matrix/'
output_dir = '/g/data/zk16/qing/tuomas/shNupr1_scRNA/02scrublet/'
counts_matrix = scipy.io.mmread(input_dir +   rundate + '/' + folder_dir + '/matrix.mtx.gz').T.tocsc()
counts_matrix = counts_matrix[:,:-10]
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

def load_genes(filename, delimiter='\t', column=0, skip_rows=0):
    gene_list = []
    gene_dict = {}

    with gzip.open(filename, 'rb') as f:
        for iL in range(skip_rows):
            f.readline()
        for l in f:
            gene = l.decode().strip('\n').split(delimiter)[column]
            if gene in gene_dict:
                gene_dict[gene] += 1
                gene_list.append(gene + '__' + str(gene_dict[gene]))
                if gene_dict[gene] == 2:
                    i = gene_list.index(gene)
                    gene_list[i] = gene + '__1'
            else: 
                gene_dict[gene] = 1
                gene_list.append(gene)
    return gene_list
    
genes = np.array(load_genes(input_dir + rundate +'/'+  folder_dir + 'features.tsv.gz', delimiter='\t', column=1))
genes = genes[:-10]
print('Number of genes in gene list: {}'.format(len(genes)))

#Initialize Scrublet object
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

#Run the default pipeline
#Doublet simulation
#Normalization, gene filtering, rescaling, PCA
#Doublet score calculation
#Doublet score threshold detection and doublet calling
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

#change threshold
scrub.call_doublets(threshold=0.25)
                                                    
scrub.plot_histogram()

scrub.plot_histogram()[0].savefig(output_dir + rundate + '_doublet_scores.pdf')

print(doublet_scores[0:5])
print(doublet_scores.shape)
pd.DataFrame(doublet_scores).to_csv(output_dir +rundate+ '_doublet_scores.txt',sep='\t')

#doublet_scores.savetxt(input_dir + '/doublet_scores.txt', doublet_scores, delimiter='\t')
#doublet_scores_file = open(input_dir + '/doublet_scores.txt', 'w')
#doublet_scores_file.write(str(doublet_scores))
#doublet_scores_file.close()

#Get 2-D embedding to visualize the results
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')

#Plot doublet predictions on 2-D embedding
scrub.plot_embedding('UMAP', order_points=True)[0].savefig(output_dir + rundate+ '_doublet_UMAP.pdf')
