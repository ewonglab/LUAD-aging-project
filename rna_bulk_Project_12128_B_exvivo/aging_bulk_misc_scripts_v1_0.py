################################################################################
######################### 10X -- MISCELLANEOUS SCRIPTS #########################
################################################################################

"""
A variety of smaller scripts for data input, data wrangling and transformation,
data plotting and data analysis.
Version: Python 3
"""

################################################################################
################################ DEPENDENCIES ##################################
################################################################################

import os, math, datetime, random, itertools
from collections import Counter
import numpy as np
import scipy.stats
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch
from fastcluster import linkage
from polo import optimal_leaf_ordering
import anndata

from aging_bulk_plot_v1_0 import *


################################################################################
################################# DATA INPUT ###################################
################################################################################

def create_ID():
    
    "Creates experiment ID (YmdHm) to identify output"
    
    exp_id = datetime.datetime.now().strftime("%Y%m%d%H%M")
    
    print("\nThe experiment ID is %s" % (exp_id))
    
    return exp_id
    
################################################################################

def save_to_txt(dataset, path, id_, name):
    
    """
    Saves pd.DataFrames or pd.Series to csv.
    ----------
    dataset: [pd.DataFrame] or [pd.Series].
    path: path to saving location.
    id_: experimental ID (e.g. YYMMDDHHMM).
    name: name of saved file. Format: /path/ID_name.
    """
            
    dataset.to_csv('%s/%s_%s.txt' % (path, id_, name), sep = '\t')
    
################################################################################

def save_to_pickle(dataset, path, id_, name):
    
    """
    Saves pd.DataFrames or pd.Series to pickle.
    ----------
    dataset: [pd.DataFrame] or [pd.Series].
    path: path to saving location.
    id_: experimental ID (e.g. YYMMDDHHMM).
    name: name of saved file. Format: /path/ID_name.
    """

        
    dataset.to_pickle('%s/%s_%s.txt' % (path, id_, name))

################################################################################

def save_to_hdf(df, path, id_, name):

    df.to_hdf('%s/%s_%s.h5' % (path, id_, name),
              key = name,
              mode = 'w',
              format = 'fixed')

################################################################################

def load_from_txt(path, id_, name, dform):
    
    """
    loads pd.DataFrames or pd.Series from csv.
    ----------
    path: path to saving location.
    id_: experimental ID (e.g. YYMMDDHHMM).
    name: name of saved file. Format: /path/ID_name.
    datatype: 'DataFrame' or 'Series'.
    ----------
    returns [pd.DataFrame] or [pd.Series]
    """
    
    if dform == 'DataFrame':
        return pd.read_csv('%s/%s_%s.txt' % (path, id_, name), sep = '\t', header = 0, index_col = 0, 
                           low_memory = False, squeeze = True)
    
    elif dform == 'Series':
        return pd.read_csv('%s/%s_%s.txt' % (path, id_, name), sep = '\t', header = None, index_col = 0, 
                           low_memory = False, squeeze = True)
    
################################################################################
    
def load_from_pickle(path, id_, name):
    
    """
    loads pd.DataFrames or pd.Series from csv.
    ----------
    path: path to saving location.
    id_: experimental ID (e.g. YYMMDDHHMM).
    name: name of saved file. Format: /path/ID_name.
    ----------
    returns [pd.DataFrame] or [pd.Series]
    """
    
    return pd.read_pickle('%s/%s_%s.txt' % (path, id_, name))

################################################################################

def load_from_hdf(path, id_, name):
    
    return pd.read_hdf('%s/%s_%s.h5' % (path, id_, name), key = name)


################################################################################
################# DATA TRANSFORMATION AND FEATURE SELECTION ####################
################################################################################

def filter_nonexpressed(df):
    
    print('%s genes in dataset' % len(df.index))
    
    genes_sel = df.sum(axis=1)[df.sum(axis=1)>0].index
    print('After removing non-expressed genes, %s genes remain' % (len(genes_sel)))

    return df.loc[genes_sel]

################################################################################

def filter_genes(df, min_sum=None, min_mean=None, min_cells=None):
    
    genes_sel = df.index
    print('%s genes in dataset' % len(genes_sel))
    
    if min_sum:
        genes_sel = df.loc[genes_sel].sum(axis=1)[df.loc[genes_sel].sum(axis=1)>=min_sum].index
        print('After removing genes with total expression of less than %s reads, %s genes remain' % (min_sum, len(genes_sel)))
        
    if min_mean:
        genes_sel = df.loc[genes_sel].mean(axis=1)[df.loc[genes_sel].mean(axis=1)>=min_mean].index
        print('After removing genes with mean expression of less than %s reads, %s genes remain' % (min_mean, len(genes_sel)))
        
    if min_cells:
        genes_sel = (df.loc[genes_sel]>0).sum(axis=1)[(df.loc[genes_sel]>0).sum(axis=1)>=min_cells].index
        print('After removing genes expressed in less than %s cells, %s genes remain' % (min_cells, len(genes_sel)))
        
    return df.loc[genes_sel]


################################################################################

def select_features_polyfit_v2(data, cutoff_mean, n_features, return_all=False):
    
    ####################
    
    def log2_var_polyfit(dataset):
    
        """
        """
        
        data_mean = dataset.mean(axis = 1)
        data_var = dataset.var(axis = 1)
        
        log2_mean = np.log2(data_mean + 1)
        log2_var = np.log2(data_var + 1)
            
        z = np.polyfit(log2_mean, log2_var, 2)
        
        log2_var_fit = z[0] * (log2_mean**2) + z[1] * log2_mean + z[2]
        log2_var_diff = log2_var - log2_var_fit
    
        return log2_var_fit, log2_var_diff, z

    ####################

    data = filter_genes(data, min_mean=cutoff_mean)
        
    log2_var_fit, log2_var_diff, z = log2_var_polyfit(data)
    
    genes_sel = log2_var_diff.sort_values()[-n_features:].index
    
    draw_log2_var_polyfit(data, log2_var_diff, z, selected=genes_sel)
    
    print("\nAfter high variance feature selection, %s genes remain" % (len(genes_sel)))
    
    data = data.loc[genes_sel]
    
    data_log2 = np.log2(data + 1 )
    
    if return_all==True:
        return data_log2, log2_var_diff, z
    
    else:
        return data_log2
    
################################################################################

def select_features_vst(df, n_features, loess_span=0.3, clip='auto'):
    
    """
    Similar to the Seurat approach found under FindVariableFeatures.
    """
    
    ################################
    
    def get_stdvar(df,mu,sd_exp,clip):
    
        df = (df - mu) / sd_exp #standardize values
        df = df.clip(upper=clip) #clip max
        df = np.square(df) #extract standardized feature variance
        return df.sum() / (len(df)-1)

    ################################
    
    from rpy2 import robjects
    from rpy2.robjects.packages import importr
    rbase = importr('base')
    rstats = importr('stats')
    
    #get mean ad
    
    hvf = pd.DataFrame(index=df.index)
    hvf['mu'] = df.mean(axis=1)
    hvf['var'] = df.var(axis=1)
    
    #exclude constant rows
    
    i = hvf['var'][hvf['var']>0].index
    df = df.loc[i]
    hvf = hvf.loc[i]
    
    #use local polynomial regression (loess, R implementation)
    
    rmu = robjects.FloatVector(hvf['mu'])
    rvar = robjects.FloatVector(hvf['var'])
    
    fit = rstats.loess(formula = robjects.Formula('log10(x = var) ~ log10(x = mu)'), 
                   data = robjects.DataFrame({'mu':rmu,'var':rvar}), 
                   span = loess_span)
    hvf['var_exp'] = [10**i for i in fit.rx2['fitted']]
    hvf['std_exp'] = np.sqrt(hvf['var_exp'])
        
    #standardize matrix based on expected standard deviation
    
    if clip == 'auto': clip = np.sqrt(len(df.columns))
        
    hvf['var_std'] = [get_stdvar(df.loc[i], hvf.loc[i,'mu'], hvf.loc[i, 'std_exp'], clip) for i in hvf.index]
    
    #return features with highest standardized variance
    
    return hvf['var_std'].sort_values(ascending = False)[:n_features].index

################################################################################

def pca_explained_var(data, dim = 50, **kwargs):
    
    pca = PCA(n_components=dim, **kwargs)
    pca_fit = pca.fit(data.T)
    exp_var = pca_fit.explained_variance_
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(range(dim), exp_var)
    
    plt.figure(figsize = [10,5], facecolor = 'w')
    ax = plt.axes()
    
    ax.set_xlabel('Principal components')
    ax.set_ylabel('Explained Variance')
    
    ax.set_xlim(-0.5, dim-0.5)
    ax.set_ylim(0, np.max(exp_var) * 1.1)
    
    ax.scatter(range(dim), exp_var, c = 'dodgerblue', linewidth = 0, s = 75)

################################################################################

def dim_reduc(df, dim=50, method='PCA',**kwargs):
    
    from sklearn.decomposition import PCA, TruncatedSVD, NMF
    
    if method == 'PCA':
        pca = PCA(n_components=dim, **kwargs)
        return pd.DataFrame(pca.fit_transform(df.T).T, index = range(dim), columns = df.columns)
    
    elif method == 'TruncatedSVD':
        tSVD = TruncatedSVD(n_components=dim, **kwargs)
        return pd.DataFrame(tSVD.fit_transform(df.T).T, index = range(dim), columns = df.columns)
    
    elif method == 'NMF':
        nmf = NMF(n_components=dim, **kwargs)
        return pd.DataFrame(nmf.fit_transform(df.T).T, index = range(dim), columns = df.columns)

################################################################################

def get_size_factors(counts):
    return np.exp2(np.log2(counts+1).apply(lambda x: x - np.mean(x), axis=1).median(axis=0))

################################################################################

def groups_reorder(groups, order, link_to = None):
    
    """
    Reorders the groups in an sample or gene group Series either completely or partially
    -----
    groups: pd.Series of either samples (Cell ID) or gene (gene ID) linked to groups (int)
    order: list containing either complete or partial new order of groups
    link_to: defines which group position is retained when groups are reorded partially; default == None, groups are linked to
    first group in 'order'
    -----
    returns reordered group Series
    """
    
    # (1) Define new group order
    
    if set(order) == set(groups):
        order_new = order
        
    else:
        
        order_new = return_unique(groups, drop_zero = False)
        
        if link_to in order:
            link = link_to
        
        elif link_to not in order or link_to == None:
            link = order[0]
            
        order.remove(link)
        
        for group in order:
            
            order_new.remove(group)
            ins_ix = order_new.index(link) + 1
            gr_ix = order.index(group)
            order_new.insert(ins_ix + gr_ix, group)
            
    # (2) Reorder groups
    
    groups_new = pd.Series()
    
    for group in order_new:
        
        groups_new = groups_new.append(groups[groups == group])
        
    groups_new = groups_new
    
    return groups_new
    
################################################################################
############################## HELPER FUNCTIONS ################################
################################################################################

def chunks(l, n):
    """ 
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]
