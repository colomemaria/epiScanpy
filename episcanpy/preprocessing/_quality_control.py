import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad
from scipy.sparse import issparse

def cal_var(adata, show=True):
    
    if issparse(adata.X):
        adata.var['n_cells'] = adata.X.sum(axis=0).tolist()[0] ## double check this
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = abs(adata.var['prop_shared_cells']-0.5)
    else:
        adata.var['n_cells'] = adata.X.sum(axis=0)
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = abs(adata.var['prop_shared_cells']-0.5)
    
    if show: # plotting
        
        fig = plt.figure(figsize=(8,6))
        sns.distplot(adata.var['variability_score'], bins=40)
        sns.distplot(adata.var['prop_shared_cells'], bins=40)
        fig.legend(labels=['variability_score','prop_shared_cells'],
                   loc='upper right')
        plt.ylabel('nb of features')
        plt.title('Distribution of feature coverage')
        plt.show()
        
def select_var_feature(adata, min_score=0.5, nb_features=None, show=True, copy=False):
    """
    This function computes a score to rank the most variable features across all cells. 
    Then it selects the most variable features according to either a specified number of features (nb_features) or a minimum variance score (min_score).
    
    adata: adata object
    
    min_score: minimum threshold of the variability score to retain features, where ### is the score of the most variable features and #### is the score of the least variable features.
    
    nb_features: default value is None, if specify it will select a the top most variable features.
    if the nb_features is larger than the total number of feature, it filters based on the min_score argument
    
    show: default value True, it will plot the distribution of var.
    
    copy: return a new adata object if copy == True.
    """
    if copy:
        inplace=False
    else:
        inplace=True

    adata = adata.copy() if not inplace else adata
    
    # calculate variability score
    cal_var(adata, show=show) # adds variability score for each feature 
    # adata.var['variablility_score'] = abs(adata.var['prop_shared_cells']-0.5)
    var_annot = adata.var.sort_values(ascending=True, by ='variability_score')

    # calculate the min score to get a specific number of feature        
    if nb_features != None and nb_features < len(adata.var_names):  
            min_score = var_annot['variability_score'][nb_features]    
     
    
    adata_tmp = adata[:,adata.var['variability_score']<=min_score].copy()
    
    ## return the filtered AnnData objet.
    if not inplace:
        adata_tmp = adata[:,adata.var['variability_score']<=min_score]
        return(adata_tmp)
    else:
        adata._inplace_subset_var(adata.var['variability_score']<=min_score)
        
def binarize(adata, copy=False):
    """convert into a binary matrix"""

    if copy:
        adata2 = adata.copy()
        adata2.X[adata2.X != 0] = 1
        return(adata2)
    else:
        adata.X[adata.X != 0] = 1
        
        
import warnings
from warnings import warn

def coverage_cells(adata, bins=50, key_added=None, log=False, binary=None, xlabel=None,
    ylabel=None, title=None, color=None, edgecolor=None, save=None):
    """
    Histogram of the number of open features (in the case of ATAC-seq data) per cell.
    """
    if key_added == None:
        key_added='nb_features'
    
    # calculate the number of features per cell
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        sum_peaks = np.sum(adata.X, axis=1).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        sum_peaks = np.sum(tmp_array.X, axis=1).tolist()
  
    
    if issparse(adata.X):
        sum_peaks = [element[0] for element in sum_peaks]        

    adata.obs[key_added] = sum_peaks
    

    
    # plotting parameters
    if xlabel ==None:
        plt.xlabel('number of features')
    else:
        plt.xlabel(xlabel)
        
    if ylabel ==None:
        plt.ylabel('number of cells')
    else:
        plt.ylabel(ylabel)
    
    if title !=None:
        plt.title(title)
        
    if color == None:   
        color='c'
    if edgecolor == None:
        edgecolor='k'

    if log:
        if 0 in sum_peaks:
            warnings.warn("""Some cells do not contain any open feature. Use epi.pp.filter_cells(adata, min_features=1) to remove these cells.""")
        
        plt.xlabel('number of features (log scale)')
        fig = plt.hist(np.log(sum_peaks), bins, color=color, edgecolor=edgecolor)
    else:
        fig = plt.hist(sum_peaks, bins, color=color, edgecolor=edgecolor)
    
    #fig = plot.get_figure()
    if save!= None:
        plt.savefig(save)
plt.show()


def commonness_features(adata, threshold=None, log=False, binary=None, key_added=None, 
    bw=0.5,xlabel=None, title=None, color=None, edgecolor=None, save=None):
    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.
    """
    
    if key_added == None:
        key_added='commonness'
        
    if xlabel ==None:
        plt.xlabel('cells sharing a feature')
    else:
        plt.xlabel(xlabel)    

    
    # calculate for each feature in how many cell it is open
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        common = np.sum(adata.X, axis=0).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        common = np.sum(tmp_array.X, axis=0).tolist()

    if issparse(adata.X):
        common = common[0]        

    # save key
    adata.var[key_added] = common

    if threshold != None:
        plt.axvline(x=threshold, color='r')
    bw_param = bw
    sns.set_style('whitegrid')
    #sns.kdeplot(np.array(adata.var[key_added]), bw=bw_param)
    
    if color == None:   
        color='c'
    if edgecolor == None:
        edgecolor='k'

    if log:
        if 0 in common:
            warnings.warn("""Some features are not open in any cell. Use epi.pp.filter_features(adata, min_cells=1) to remove these features.""")
        
        plt.xlabel('cells sharing a feature (log scale)')
        fig = plt.hist(np.log(common), bins=int(100), color=color, edgecolor=edgecolor)
    else:
        fig = plt.hist(common, bins=int(80), color=color, edgecolor=edgecolor)
        
    
    if title !=None:
        plt.title(title)
    
    #fig = plot.get_figure()
    if save!= None:
        fig.savefig(save)
    plt.show()

    adata.var[key_added] = common