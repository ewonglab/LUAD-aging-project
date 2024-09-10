################################################################################
############################## 10X -- PLOT #####################################
################################################################################

"""
Scripts for drawing.
Version: Python 3
"""

################################################################################
################################ DEPENDENCIES ##################################
################################################################################

import random
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter
from mpl_toolkits.axes_grid1 import make_axes_locatable

################################################################################
############################### MISC FUNCTIONS #################################
################################################################################

def remove_ticks(axes, linewidth = 0.5):
    """
    Removes ticks from matplotlib Axes instance
    """
    axes.set_xticklabels([]), axes.set_yticklabels([])
    axes.set_xticks([]), axes.set_yticks([])
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(linewidth)
        
################################################################################


def return_unique(groups, drop_zero = False):
    """
    Returns unique instances from a list (e.g. an AP cluster Series) in order 
    of appearance.
    """
    unique = []
    
    for element in groups.values:
        if element not in unique:
            unique.append(element)
            
    if drop_zero == True:
        unique.remove(0)
        
    return unique

################################################################################

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


################################################################################
################################## PLOTTING ####################################
################################################################################

def draw_polyfit(dataset, log2_var_diff, z, selected):
    
    """
    """
    
    data_mean = dataset.mean(axis = 1)
    data_var = dataset.var(axis = 1) 
    
    log2_mean = np.log2(data_mean + 1)
    log2_var = np.log2(data_var + 1)
    
    line_x = np.arange(log2_mean.min(), log2_mean.max(), 0.01)
    line_y = [z[0] * (x**2) + z[1] * x + z[2] for x in line_x]
    
    clist = pd.Series('blue', index = log2_mean.index)
    clist[log2_var_diff[log2_var_diff > 0].index] = 'red'
    
    if np.all(selected != None):
        clist[selected] = 'green'
        
    plt.figure(figsize = [10,10], facecolor = 'w')
    ax0 = plt.axes()
    
    ax0.set_xlabel('Mean number of transcripts [log2]')
    ax0.set_ylabel('Variance [log2]')
    
    ax0.set_xlim(log2_mean.min() - 0.5, log2_mean.max() + 0.5)
    ax0.set_ylim(log2_var.min() - 0.5, log2_var.max() + 0.5)
    
    ax0.scatter(log2_mean, log2_var, c = clist, linewidth = 0,)
    ax0.plot(line_x, line_y, c = 'r', linewidth = 3)
    
################################################################################

def draw_heatmap(df, cols, rows, cmap_c, cmap_r):

    df = df.ix[list(rows.index), list(cols.index)]
    df = df.apply(lambda x: x / max(x), axis = 1)

    plt.figure(figsize=(20,20), facecolor = 'w')
    
    #draw heatmap

    ax0 = plt.axes()
    ax0.set_position([0.05, 0.05, 0.9, 0.9])
    
    ax0.imshow(df, aspect = 'auto', interpolation = 'nearest')
    
    remove_ticks(ax0)
    
    #draw row bar

    divider = make_axes_locatable(ax0)
    ax2 = divider.append_axes("right", size= 0.5, pad=0.05)

    ax2.set_xlim(0,1)
    ax2.set_ylim(len(rows),0)
    
    for p, i in enumerate(rows):
        ax2.axhspan(p, p+1, color = cmap_r[i])

    remove_ticks(ax2)
    
    #draw row ticks
    
    row_ticks = pd.Series(index = set(rows))
    
    for gr in row_ticks.index:
                
        first_ix = list(rows.values).index(gr)
        last_ix = len(rows) - list(rows.values)[::-1].index(gr)
        row_ticks[gr] = first_ix + ((last_ix - first_ix) / 2.0)
        
    ax2.set_yticks(row_ticks.values)
    ax2.set_yticklabels(row_ticks.index)
    ax2.yaxis.set_ticks_position('right')
    
    #draw col bar
    
    ax3 = divider.append_axes("bottom", size= 0.5, pad=0.05)

    ax3.set_ylim(0, 1)
    ax3.set_xlim(0, len(cols))
    
    for p, i in enumerate(cols):
        ax3.axvspan(p, p+1, color = cmap_c[i])
    
    remove_ticks(ax3)
    
    #draw col ticks
    
    col_ticks = pd.Series(index = set(cols))
        
    for gr in col_ticks.index:
                
        first_ix = list(cols.values).index(gr)
        last_ix = len(cols) - list(cols.values)[::-1].index(gr)
        col_ticks[gr] = first_ix + ((last_ix - first_ix) / 2.0)
        
    ax3.set_xticks(col_ticks.values)
    ax3.set_xticklabels(col_ticks.index)
    ax3.xaxis.set_ticks_position('bottom')

################################################################################

def draw_dist_mat(dist_mat, groups, **kwargs):
    
    """
    Draws distance matrices of either m cells or n genes randomly shuffled and
    ordered according to group Series (e.g. AP clustering).
    ----------
    dist_mat: pd.DataFrame with distances of either m x m cells or n x n genes.
    groups: pd.Series with ordered cluster identity of m cells or n genes.
    """
    
    plt.figure(figsize = [20,10], facecolor = 'w')

    ax0 = plt.subplot(121)

    tmp_ix = list(dist_mat.index)
    random.shuffle(tmp_ix)

    ax0.matshow(dist_mat.ix[tmp_ix, tmp_ix], **kwargs)

    ax1 = plt.subplot(122)

    ax1.matshow(dist_mat.ix[groups.index, groups.index], **kwargs)
    
################################################################################
    
def draw_barplots(dataset, cell_groups, genes, cmap = plt.cm.jet):
    
    """
    draws expression of selected genes in order barplot with juxtaposed group identity
    dataset: pd.DataFrame of n samples over m genes
    sample_group_labels: ordered (!) pd.Series showing sample specific group indentity 
    list_of_genes: list of selected genes
    color: matplotlib cmap
    """
    
    # set figure framework
    
    plt.figure(facecolor = 'w', figsize = (21, len(genes) * 3 + 1))
        
    gs0 = plt.GridSpec(1,1, left = 0.05, right = 0.95, top = 1 - 0.05 / len(genes),
                       bottom = 1 - 0.15 / len(genes), hspace = 0.0, wspace = 0.0)
    
    gs1 = plt.GridSpec(len(genes), 1, hspace = 0.05, wspace = 0.0, left = 0.05, right = 0.95, 
                       top = 1 - 0.2 / len(genes) , bottom = 0.05)
    
    #define dataset
    
    dataset = dataset.loc[genes, cell_groups.index]

     #define colormap
    
    if type(cmap) != dict:
        cm = cmap
        cmap = {}
        for ix, gr in enumerate(return_unique(cell_groups)):
            cmap[gr] = cm(float(ix) / len(set(cell_groups)))
            
    clist = [cmap[cell_groups[c]] for c in cell_groups.index]
    
    #draw genes
    
    for ix, g in enumerate(genes):
    
        ax = plt.subplot(gs1[ix])
        ax.set_xlim(left = 0, right = (len(dataset.columns)))
                     
        if g != genes[-1]:
            ax.xaxis.set_ticklabels([])
        
        elif g == genes[-1]:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(15) 
                
        ax.set_ylabel(g, fontsize = 25)
        ax.yaxis.labelpad = 10
        
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(15)
            tick.label.set_label(10)
            
        ax.bar(np.arange(0, len(dataset.columns),1), 
               dataset.loc[g],
               width=1,
               color=clist,
               linewidth=0)
    
    #create groups bar
    
    ax_bar = plt.subplot(gs0[0])
    
    ax_bar.set_xlim(left = 0, right = (len(dataset.columns)))
    ax_bar.set_ylim(bottom = 0, top = 1)
    
    for ix, val in enumerate(cell_groups.values):
        
        ax_bar.axvspan(ix,
                       ix+1,
                       ymin=0,
                       ymax=1, 
                       color = clist[ix])
        
    remove_ticks(ax_bar)
    
################################################################################
    
def QC_hist(data, param, batches, cmap_batches = plt.cm.tab20, xlim = (0,25000), **kwargs):
    
    #compile data based on param
    
    if param == 'reads': data = data.sum(axis=0)
    elif param == 'genes': data = (data>0).sum(axis=0)
            
    #initialize figure
    
    l = len(set(batches))

    height = 2 * l
    width = 7.5
    plt.figure(facecolor = 'w', figsize = (width, height))
    gs = plt.GridSpec(l, 1, hspace=0)
    
    #get colormap
    
    if type(cmap_batches) != dict:
        cm = cmap_batches
        cmap = {}
        for i, b in enumerate(set(batches)):
            cmap[b] = cm(float(i) / 20)
            
    else:
        cmap = cmap_batches
    
    #plot batches
    
    for i, b in enumerate(return_unique(batches)):
        
        ax = plt.subplot(gs[i])
        
        ax.set_xlim(xlim)
        
        if i == l-1:
            if param == 'reads': ax.set_xlabel('Number of reads per cell')
            elif param == 'genes': ax.set_xlabel('Number of genes per cell')
                
        else:
            ax.set_xticks([])
                
        ax.set_ylabel('Number of cells')
        
        i_tmp = batches[batches==b].index
        
        ax.hist(data[i_tmp],
            bins=100,
            color=cmap[b],
            range = xlim,
            **kwargs)
            
        ax.axvline(np.median(data[i_tmp]), color = 'black', linestyle = '-', linewidth = 3)
        ax.axvline(np.percentile(data[i_tmp], 5), color = 'black', linestyle = '--', linewidth = 2)
        ax.axvline(np.percentile(data[i_tmp], 95), color = 'black', linestyle = '--', linewidth = 2)
        ax.axvline(np.percentile(data[i_tmp], 1), color = 'black', linestyle = ':', linewidth = 2)
        ax.axvline(np.percentile(data[i_tmp], 99), color = 'black', linestyle = ':', linewidth = 2)
        
        ax.text(ax.get_xlim()[1] * 0.95, ax.get_ylim()[1] * 0.8, b, fontsize = 15, ha = 'right')

################################################################################

def QC_bar(data,param,meta,cmap,**kwargs):
    
    #compile data based on param
    
    if param == 'reads': data = data.sum(axis=0).sort_values()
    elif param == 'genes': data = (data>0).sum(axis=0).sort_values()
        
    #plot data
    
    #initialize figure

    height = 5
    width = 10
    plt.figure(facecolor = 'w', figsize = (width, height))
    
    ax = plt.subplot(111)
    ax.set_xlim(0, len(data.index))
    if param == 'reads': ax.set_ylabel('Number of reads per cell')
    elif param == 'genes': ax.set_ylabel('Number of genes per cell')
        
    #defined colors
    
    clist = [cmap[meta[i]] for i in data.index]
    
    #plot bars
    
    ax.bar(range(len(data.index)),
           data, color = clist,
           width = 1.0,
           **kwargs)
    
    #plot median and mean
    
    ax.axhline(np.median(data), color = 'black', linestyle = '-', linewidth = 3)
    ax.axhline(np.percentile(data, 5), color = 'black', linestyle = '--', linewidth = 2)
    ax.axhline(np.percentile(data, 95), color = 'black', linestyle = '--', linewidth = 2)
    ax.axhline(np.percentile(data, 1), color = 'black', linestyle = ':', linewidth = 2)
    ax.axhline(np.percentile(data, 99), color = 'black', linestyle = ':', linewidth = 2)

################################################################################

def draw_scatter_groups(coords, groups, cmap = plt.cm.tab20, pad = 2, s = 50, 
                        show_axes = True, show_legend=True):
    
    #initialize figure

    height = 7
    width = 10

    plt.figure(facecolor = 'w', figsize = (width, height))
    gs = plt.GridSpec(1,2, wspace=0.025, width_ratios=[7,3])

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
        xlim = (x_cent - 0.5 * y_diff - pad, x_cent + 0.5 * y_diff + pad,)
        ylim = (y_min - pad, y_max + pad)

    #define x- and y-axes

    ax = plt.subplot(gs[0])

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    
    #define colormap
    
    if type(cmap) != dict:
        cm = cmap
        cmap = {}
        for ix, gr in enumerate(return_unique(groups)):
            cmap[gr] = cm(float(ix) / 20)
            
    clist = [cmap[groups[c]] for c in groups.index]
    
    #plot

    ax.scatter(coords[:,0],
               coords[:,1], 
               s = s,
               linewidth = 0.0,
               c = clist)
    
    if not show_axes:
        clean_axis(ax)
    
    #plot legend
    
    if show_legend:
    
        grs = list(set(groups))
        grs.sort()

        ax = plt.subplot(gs[1])
        ax.set_xlim(0,1)
        if len(grs) <= 10: ax.set_ylim(10.5, -0.5)
        else: ax.set_ylim(len(grs) + 0.5, -0.5)

        for p, i in enumerate(grs):
            ax.scatter(0.15, p, color = cmap[i], s = 200)
            ax.text(0.3, p, i, fontsize = 15, family = 'Arial', va = 'center')

        clean_axis(ax)

################################################################################

def draw_scatter_expr(coords, expr, vmin, vmax,text = None,  cmap = plt.cm.viridis, pad = 2, s = 50, 
                      show_axes = True, show_legend=True):
    
    if vmax < 1: 
        vmax = 1
    
    #initialize figure
    
    height = 7
    width = 7.5
    plt.figure(facecolor = 'w', figsize = (width, height))
    gs = plt.GridSpec(1,2, wspace=0.025, width_ratios=[7,.5])

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
        xlim = (x_cent - 0.5 * y_diff - pad, x_cent + 0.5 * y_diff + pad,)
        ylim = (y_min - pad, y_max + pad)

    #define x- and y-axes

    ax = plt.subplot(gs[0])

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    
    #define colormap
            
    clist = np.array([cmap((e-vmin)/(vmax-vmin)) for e in expr])
    ixz = expr.argsort()
        
    #plot

    ax.scatter(coords[:,0][ixz],
               coords[:,1][ixz], 
               s = s,
               linewidth = 0.0,
               c = clist[ixz.values])
    
    #text
    
    if text:
        ax.text(xlim[0] + (xlim[1]-xlim[0]) / 2,
                ylim[1] * 1.05,
                text,
                fontsize = 15, va = 'center', ha = 'center')
        
    if not show_axes:
        clean_axis(ax)
    
    #plot colorbar
    
    if show_legend:

        ax = plt.subplot(gs[1])

        ax.set_xlim(0,1)
        ax.set_xticks([])

        ax.set_ylim(vmin, vmax)
        ax.yaxis.set_ticks_position('right')

        for i in np.linspace(vmin, vmax, 100):
            ax.axhspan(i, i + (vmax-vmin) / 100, color = cmap((i-vmin)/(vmax-vmin)))

################################################################################

def draw_fraction_cluster(cluster, variable, cmap, cl_order=False, var_order=False):
    
    if not cl_order:
        cl_order = return_unique(cluster)
        
    if not var_order:
        var_order = return_unique(variable)
        
    #initialize figure

    height = 10
    width = len(set(cluster)) * 1.5
    
    fig = plt.figure(facecolor = 'w', figsize = (width, height))

    #set axes

    ax = plt.subplot(111)

    ax.set_xlim(-0.5, len(cl_order)-0.5)
    ax.set_xticks(list(range(len(cl_order))))
    ax.set_xticklabels(cl_order, fontsize = 30, family = 'Arial', rotation = 'vertical')

    ax.set_ylim(0,1)
    ax.set_yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    ax.set_yticklabels([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], family = 'Arial', fontsize = 25)
    ax.set_ylabel('Fraction', fontsize = 40, family = 'Arial')
    
    #iterate over clusters
    
    for pos, cl in enumerate(cl_order):
        c_tmp = cluster[cluster==cl].index
        l = len(c_tmp)
        rep_tmp = Counter(variable[c_tmp])
        y = 0
        for var in var_order:
            y_new = y + rep_tmp[var]/l
            ax.bar(x = pos, bottom=y, height=y_new, width=1, color = cmap[var])
            y = y_new
        ax.axvline(pos+0.5, color = 'k', linewidth = 0.5)

################################################################################

def draw_volcano_plot_diffxpy(de_sum, highlights = None, xlim = None, ylim = None, sign_thr = -np.log10(1e-10)):
    
    "Use summary document returned by diffxpy"
    
    #initialize figure

    height = 5
    width = 5

    plt.figure(facecolor = 'w', figsize = (width, height))
    
    fc = de_sum['log2fc']
    pval = -np.log10(de_sum['qval'])
        
    #set infinite pval to max observed non-infinite pval 

    pval_max = np.max(np.ma.array(pval, mask=np.isinf(pval))) #max on masked array
    pval[np.isinf(pval)] = pval_max
        
    #initialize axis and set axes limits
    
    ax = plt.subplot(111)
    
    ax.set_xlabel('log2(Fold change)', family = 'Arial', fontsize = 20)
    ax.set_ylabel('-log10(P-value)', family = 'Arial', fontsize = 20)
    
    if xlim: ax.set_xlim(xlim)
    else: ax.set_xlim(fc.min() * 1.1, fc.max() * 1.1)
        
    for tick_pos, tick in enumerate(ax.xaxis.get_major_ticks()):
        ax.xaxis.get_major_ticks()[tick_pos].label.set_family('Arial')
        ax.xaxis.get_major_ticks()[tick_pos].label.set_fontsize(15)
        
    if ylim: ax.ylim(ylim)
    else: ax.set_ylim(0, pval_max * 1.1)
        
    for tick_pos, tick in enumerate(ax.yaxis.get_major_ticks()):
        ax.yaxis.get_major_ticks()[tick_pos].label.set_family('Arial')
        ax.yaxis.get_major_ticks()[tick_pos].label.set_fontsize(15)
                
    #plot data
    
    ax.scatter(fc, pval, s = 10, color = 'silver')
    
    #plot highlighted genes
    
    if highlights:
        ax.scatter(fc[highlights], pval[highlights], s = 50, color = 'red')
    
    #plot significance threshold
    
    ax.axhline(sign_thr, linewidth = 1, color = 'red')

################################################################################