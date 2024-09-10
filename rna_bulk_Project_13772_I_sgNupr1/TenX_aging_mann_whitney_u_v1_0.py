################################################################################
############################ 10X -- MANN-WHITNEY U #############################
################################################################################

"""
Scripts for non-parametric statistical testing using Mann-Whitney U test.
Version: Python3
"""

################################################################################
################################ DEPENDENCIES ##################################
################################################################################

import random, itertools
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu

import rpy2.robjects as robj
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import Vector, DataFrame, FloatVector, IntVector, StrVector, ListVector, Matrix, BoolVector
#from rpy2.rinterface import RNULLType
stats = importr('stats')

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

def mwu_vs_average(data, groups, genes, dview, BH = True, log = True):
        
    #########################
    
    def mwu_vs_average_helper(data, groups, gene):
        
        output = pd.DataFrame(index = [gene], columns = return_unique(groups))
                
        for gr in return_unique(groups):
            d1 = data.loc[groups[groups==gr].index]
            d2 = data.loc[groups[groups!=gr].index]
             
            try:
                output.loc[gene,gr] = mwu(d1, d2, alternative = 'greater')[1]
            except:
                output.loc[gene,gr] = 1.0
                
        return output
            
    #########################
    
    l = len(genes)
    
    output_tmp = dview.map_sync(mwu_vs_average_helper,
                                [data.loc[g] for g in genes], 
                                [groups] * l, 
                                genes)
        
    output = pd.concat(output_tmp, axis = 0)
        
    if BH == True:
        for col in output.columns:
            output[col] = stats.p_adjust(FloatVector(output[col]), method = 'BH')
            
    if log == True:
        output = -np.log10(output.astype(float))
    
    return output
    
################################################################################

def mwu_vs_groups(data, groups, genes, dview, BH = True, log = True):
        
    #########################
    
    def mwu_vs_groups_helper(data, groups, gene):
        
        output = pd.DataFrame(index = [gene], columns = return_unique(groups))
        
        for gr1 in return_unique(groups):
            d1 = data.loc[groups[groups==gr1].index]
            pvals = []
            
            for gr2 in [gr2 for gr2 in return_unique(groups) if gr2 != gr1]:
                d2 = data.loc[groups[groups==gr2].index]
                
                try:
                    pval_tmp = mwu(d1, d2, alternative = 'greater')[1]
                except:
                    pval_tmp = 1.0
                                        
                pvals.append(pval_tmp)
                    
            output.loc[gene, gr1] = np.max(pvals)
                
        return output.astype(float)
            
    #########################
    
    l = len(genes)
    
    output_tmp = dview.map_sync(mwu_vs_groups_helper,
                                [data.loc[g] for g in genes], 
                                [groups] * l, 
                                genes)
    
    output = pd.concat(output_tmp, axis = 0)
    
    if BH == True:
        for col in output.columns:
            output[col] = stats.p_adjust(FloatVector(output[col]), method = 'BH')
            
    if log == True:
        output = -np.log10(output.astype(float))
    
    return output

################################################################################

def log2_fold_change_vs_average(data, groups, genes, dview):
        
    #########################
    
    def log2_fold_change_vs_average_helper(data, groups, gene):
        
        output = pd.DataFrame(index = pd.MultiIndex.from_tuples([gene]), 
                              columns = return_unique(groups))
                
        for gr in return_unique(groups):
            d1 = data.loc[groups[groups==gr].index].values
            d2 = data.loc[groups[groups!=gr].index].values
             
            output.loc[gene,gr] = np.mean(d1) - np.mean(d2)
                
        return output
            
    #########################
    
    l = len(genes)
    
    output_tmp = dview.map_sync(log2_fold_change_vs_average_helper,
                                [data.loc[g] for g in genes], 
                                [groups] * l, 
                                genes)
        
    output = pd.concat(output_tmp, axis = 0)
    
    return output

################################################################################

def log2_fold_change_vs_groups(data, groups, genes, dview):
        
    #########################
    
    def log2_fold_change_vs_groups_helper(data, groups, gene):
        
        output = pd.DataFrame(index = pd.MultiIndex.from_tuples([gene]),
                              columns = return_unique(groups))
        
        for gr1 in return_unique(groups):
            d1 = data.loc[groups[groups==gr1].index].values
            fcs = []
            
            for gr2 in [gr2 for gr2 in return_unique(groups) if gr2 != gr1]:
                d2 = data.loc[groups[groups==gr2].index].values
                
                fc_tmp = np.mean(d1) - np.mean(d2)                         
                fcs.append(fc_tmp)
                    
            output.loc[gene, gr1] = np.min(fcs)
                
        return output.astype(float)
            
    #########################
    
    l = len(genes)
    
    output_tmp = dview.map_sync(log2_fold_change_vs_groups_helper,
                                [data.loc[g] for g in genes], 
                                [groups] * l, 
                                genes)
    
    output = pd.concat(output_tmp, axis = 0)

    return output
    
################################################################################

def mwu_select_features(data, groups, cutoff_mean):
    
    #get mean data per group
    
    mean = pd.DataFrame(index = data.index, columns = set(groups))
    
    for gr in set(groups):
        c_sel = groups[groups==gr].index
        
        mean.loc[data.index, gr] = data.loc[data.index, c_sel].mean(axis = 1)
        
    #select genes over cutoff
    
    genes_sel = mean.max(axis=1)[mean.max(axis=1)>cutoff_mean].index
    
    return genes_sel

################################################################################