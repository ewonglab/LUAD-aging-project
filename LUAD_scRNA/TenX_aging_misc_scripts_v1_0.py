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
#from polo import optimal_leaf_ordering
import anndata

from TenX_aging_plot_v1_0 import *


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

def pd2sc(df):
    
    X = df
    
    obs = pd.DataFrame(index = X.columns)
    
    var = pd.DataFrame(index = X.index)
    
    scdata = sc.AnnData(np.array(X.T), obs = obs, var = var)
    scdata.var_names_make_unique()
    scdata.obs_names_make_unique()

    return scdata

################################################################################

def sc2pd(scdata, layer, sparse=True):

    if sparse:
    
        if layer=='X':
            return pd.DataFrame(scdata.X.T.todense(),
                                index = scdata.var_names,
                                columns = scdata.obs_names)
        else:
            return pd.DataFrame(scdata.layers[layer].T.todense(),
                                index = scdata.var_names,
                                columns = scdata.obs_names)

    else:
    
        if layer=='X':
            return pd.DataFrame(scdata.X.T,
                                index = scdata.var_names,
                                columns = scdata.obs_names)
        else:
            return pd.DataFrame(scdata.layers[layer].T,
                                index = scdata.var_names,
                                columns = scdata.obs_names)

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

def filter_cells(df, min_reads=None, min_genes=None):
    
    cells_sel = df.columns
    print('%s cells in dataset' % len(cells_sel))
    
    if min_reads:
        cells_sel = df[cells_sel].sum(axis=0)[df[cells_sel].sum(axis=0)>=min_reads].index
        print('After removing cells with less than %s reads, %s cells remain' % (min_reads, len(cells_sel)))
        
    if min_genes:
        cells_sel = (df[cells_sel]>0).sum(axis=0)[(df[cells_sel]>0).sum(axis=0)>=min_genes].index
        print('After removing cells with less than %s genes, %s cells remain' % (min_genes, len(cells_sel)))
        
    return df[cells_sel]

################################################################################

def log2_transform(dataset, add = 1):
    
    """
    Calculates the binary logarithm (log2(x + y)) for every molecule count / cell x in dataset. 
    Unless specified differently, y = 1.
    ----------
    dataset: seq dataset [pd.DataFrame] containing m cells x n genes.
    add: y [float or int] in (log2(x + y)).
    ----------
    returns seq dataset [pd.DataFrame] containing m cells x n genes.
    """
    
    print('\nCalculating binary logarithm of x + %s' % (add))
    dataset = np.log2(dataset.astype(float) + float(add))
    
    return dataset

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
        df = df.clip(upper=clip) #clip maxf
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

def remove_cells(data, genes, cutoff):
    
    c_rem = []
    
    for g, c in zip(genes, cutoff):
        c_rem += list(data.loc[g][data.loc[g]>=c].index)
        
    return list(set(c_rem))

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

def order_incluster_POLO(dist, clusters, method = 'single'):
    
    ordered = []
    
    #iterate over clusters
    
    for cl in return_unique(clusters):
        
        c_sel = clusters[clusters==cl].index
        D = pdist(dist.loc[c_sel, c_sel])
        Z = linkage(D, method = method)
        optimal_Z = optimal_leaf_ordering(Z, D)
        leaves = sch.dendrogram(optimal_Z, no_plot = True)['leaves']
        ordered += list(c_sel[leaves])
        
    return clusters[ordered]

################################################################################

def calculate_complexity(anndata, plot = True, res_thr = 0.5, inplace = True, lowess_frac = 0.5):
    
    from scipy.stats import gaussian_kde
    from statsmodels.nonparametric.smoothers_lowess import lowess
    
    assert 'log1p_total_counts' in anndata.obs_keys(), 'total_counts not found, run sc.pp.calculate_qc_metrics'
    assert 'log1p_n_genes_by_counts'in anndata.obs_keys(), 'n_genes_by_counts not found, run sc.pp.calculate_qc_metrics'
    
    y = anndata.obs['log1p_total_counts'].values
    x = anndata.obs['log1p_n_genes_by_counts'].values
    
    #smooth data using lowess
    
    lowess = lowess
    sm = lowess(y,x,return_sorted=False, frac=lowess_frac)
    
    #calculate res
    
    res = y - sm
        
    #plot data
    
    if plot:
        
        ## initialize figure
        
        height = 5
        width = 5
        plt.figure(facecolor = 'w', figsize = (width, height))
        
        ## set axes
        
        ax = plt.subplot(111)
        
        ax.set_xlim(np.amin(x), np.amax(x))
        ax.set_xlabel('Number of unique genes [log1p]')
        
        ax.set_ylim(np.amin(y), np.amax(y))
        ax.set_ylabel('Number of molecules [log1p]')
                
        ## Calculate the point density
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        
        ## plot scatter
        
        idz = z.argsort()
        ax.scatter(x[idz], y[idz], c=z[idz], s=10, edgecolor='')
        
        idx = x.argsort()
        ax.plot(x[idx],sm[idx], linewidth=2, color = 'red')
        
        ## plot residuals
        
        if res_thr:     
            assert type(res_thr) is float, 'Set residual threshold (float)'
            ax.scatter(x[res>=res_thr], y[res>=res_thr], c='red', s=10, edgecolor='')
            
    if inplace:
        anndata.obs['complexity_res'] = res
    else:
        return pd.Series(res, anndata.obs_names)

################################################################################

def pseudobulk(dat, samples, metric='sum', layer=None):
    
    assert type(dat) in (pd.DataFrame, anndata._core.anndata.AnnData), 'Incorrect input data'
    assert type(samples) == pd.Series, 'Incorrect input data'
    
    if type(dat) == anndata._core.anndata.AnnData:
        assert layer, 'Specify layer'
        dat = sc2pd(dat, layer = layer)
        
    out = pd.DataFrame(index = dat.index, columns = return_unique(samples))
    
    for s in return_unique(samples):
        cs = samples[samples==s].index
        if metric == 'sum': out[s] = dat[cs].sum(axis=1)
        elif metric == 'mean': out[s] = dat[cs].mean(axis=1)
        elif metric == 'median': out[s] = dat[cs].median(axis=1)
    
    return out

################################################################################

def assign_hashes(dat, thr, log=True):
    
    #transform matrix    
    if log: dat = np.log2(dat+1)
    
    #threshold data    
    dat_bin  = dat.apply(lambda x: x > thr, axis = 0)
    
    #create output   
    return pd.Series(['-'.join(dat_bin[i][dat_bin[i]].index) if np.sum(dat_bin[i]) else 'Negative' for i in dat.columns], 
                     index = dat.columns)

################################################################################

def calculate_cell_density(coords, groups, kernel='gaussian', bw=0.5, n_bins=500, pad=2, quant_trans=True):
    
    from sklearn.neighbors import KernelDensity
    from sklearn.preprocessing import quantile_transform
    
    #define x- and y-limits

    x_min, x_max = np.min(coords[:,0]), np.max(coords[:,0])
    y_min, y_max = np.min(coords[:,1]), np.max(coords[:,1])
    x_diff, y_diff = x_max - x_min, y_max - y_min
    x_cent, y_cent = x_min + 0.5 * x_diff, y_min + 0.5 * y_diff,

    pad = pad

    if x_diff > y_diff:
        xlim = (x_min - pad, x_max + pad)
        ylim = (y_cent - 0.5 * x_diff - pad, y_cent + 0.5 * x_diff + pad,)

    if x_diff < y_diff:
        xlim = (x_cent - 0.5 * y_diff - pad, x_cent + 0.5 * y_diff + pad)
        ylim = (y_min - pad, y_max + pad)
        
    #create grid for scoring
    
    x = np.linspace(xlim[0], xlim[1], n_bins)
    y = np.linspace(ylim[0], ylim[1], n_bins)
    X,Y = np.meshgrid(y, x)
    xy = np.vstack([Y.ravel(), X.ravel()]).T
    
    #fit KDE over complete dataset to define mask
    
    kde = KernelDensity(kernel=kernel, bandwidth=bw).fit(coords)    
    Z = np.exp(kde.score_samples(xy))
    mask = np.isclose(Z,0)
    
    #define ouput
    
    out = pd.DataFrame(index = tuple(xy), columns = return_unique(groups))
    
    #fit KDE over groups
    
    for g in return_unique(groups):
        kde = KernelDensity(kernel=kernel, bandwidth=bw).fit(coords[groups==g])    
        Z = np.exp(kde.score_samples(xy))
        Z[mask] = None #set none to not condsider in scaling
        out[g] = Z
        
    #quantile transform
    
    if quant_trans:
        out = pd.DataFrame(quantile_transform(out), index = out.index, columns = out.columns)
        
    #out = out.fillna(0)
    
    return out

################################################################################

def calculate_transcriptional_distance(dat, groups, n_iter=100, n_sample=50, replace=False, distance='euclidean'):
    
    def get_distance_helper(dat1, dat2, distance='euclidean'):
        
        if distance=='euclidean':
            return scipy.spatial.distance.euclidean(dat1, dat2)
        elif distance=='correlation':
            return scipy.spatial.distance.correlation(dat1, dat2)
    
    #calculate randomly sampled pseudobulk vectors for each group
    
    psb = []
    
    for gr in set(groups):
        c_sel = groups[groups==gr].index
        psb_gr = pd.DataFrame(index = dat.index, columns = range(n_iter))
        
        for i in range(n_iter):
            c_sel_rand = np.random.choice(c_sel, n_sample, replace=replace)
            psb_gr[i] = dat[c_sel_rand].mean(axis=1)
        psb.append(psb_gr)
        
    #calculate pairwise distances between and within groups
    
    dist = pd.DataFrame(index = ['XY','XX','YY'], columns = range(n_iter))
    
    for i in range(n_iter):
        ixs = np.random.choice(range(n_iter),2, replace=False)
        dist.loc['XY'][i] = get_distance_helper(psb[0][ixs[0]], psb[1][ixs[1]], distance=distance)
        dist.loc['XX'][i] = get_distance_helper(psb[0][ixs[0]], psb[0][ixs[1]], distance=distance)
        dist.loc['YY'][i] = get_distance_helper(psb[1][ixs[0]], psb[1][ixs[1]], distance=distance)
        
    #return distance between groups normalized by median / distance within group
    
    return dist.loc['XY'] / dist.loc[['XX','YY']].median()

################################################################################

def calculate_transcriptional_distance_subgroups(dat, groups, classes, n_iter=100, 
                                                 n_sample=50, replace=False, distance='euclidean'):
    
    #define output
    
    out = pd.DataFrame(index = return_unique(groups), columns = range(n_iter))
    
    for gr in return_unique(groups):
        c_sel = groups[groups==gr].index
        out.loc[gr] = calculate_transcriptional_distance(dat, 
                                                         classes[c_sel], 
                                                         n_iter, 
                                                         n_sample, 
                                                         replace, 
                                                         distance)
    return out

################################################################################

def calculate_composition(comp, groups, sample=False, add_one=False, return_counts=False):
    
    #define output
    
    out = pd.DataFrame(index=return_unique(groups),
                       columns=return_unique(comp))
    
    #iterate over groups
    
    for g in return_unique(groups):
        i = groups[groups==g].index
        if sample:
            i = np.random.choice(i, sample, replace=False)
        cnt = Counter(comp[i])
        out.loc[g] = pd.Series(cnt)
    
    if add_one:
        out = out.fillna(1) # to allow ILR
    else:
        out = out.fillna(0)
    
    if return_counts:
        return out
    
    #calculate fractions
    
    out = out.apply(lambda x: x/np.sum(x), axis = 1)
        
    return out

################################################################################

def lda(dat, classes, class_order, return_transformation=False, **kwargs):
    
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    
    lda = LinearDiscriminantAnalysis(**kwargs)
    dat_trans = lda.fit_transform(dat, classes)
    coef = lda.coef_
    
    if class_order[0] != lda.classes_[0]:
        coef *= -1
    
    if return_transformation:
        return dat_trans, coef
    
    else:
        return coef

################################################################################

def calculate_lda_coef_sampling(comp, groups, classes, class_order, sample, n_iter=1000):
    
    #define output
    
    out = pd.DataFrame(index=range(n_iter),
                       columns=return_unique(comp))
    
    for i in range(n_iter):
        comp_i = calculate_composition(comp, groups, sample)
        coef_i = lda(comp_i, classes, class_order)
        out.loc[i] = coef_i[0]
        
    return out

################################################################################

def sample_phenotypic_volume(dat, grps, n_it, dview, sample=1000):
    
    ####################
    
    def get_phenotypic_volume(dat):
    
        cov = np.cov(dat.T)
        eigval = scipy.linalg.eigvalsh(cov)
        phenovolume = np.sum(np.log(eigval[~np.isclose(eigval,0)]))

        return phenovolume
    
    ####################
    
    out = pd.DataFrame(index=return_unique(grps), 
                       columns=range(n_it))
    
    for g in return_unique(grps):
        c_gr = grps[grps==g].index
        
        for i in range(n_it):
            c_it = np.random.choice(c_gr, size=sample, replace=False)
        
            out.loc[g, i] = get_phenotypic_volume(dat[c_it])
        
    return out

################################################################################
############################## HELPER FUNCTIONS ################################
################################################################################

def chunks(l, n):
    """ 
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]
