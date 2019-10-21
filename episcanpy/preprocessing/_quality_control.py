import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

def commonness_features(adata, threshold=None, bw=0.5,
    key_added=None, xlabel=None, title=None, color=None,
    edgecolor=None, save=None, log=False):
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
    
    common = np.sum(adata.X, axis=0).tolist()
    if len(common) == 1:
        common = [item for sublist in common for item in sublist]
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

    


def coverage_cells(adata, bins=50, key_added=None, xlabel=None,
    ylabel=None, title=None, log=False, color=None,
    edgecolor=None, save=None):
    """
    Histogram of the number of features with an open peak (in the case of ATAC-seq data)
    for every cells.
    """
    if key_added == None:
        key_added='nb_features'
    # make sum peaks
    
    #sum_peaks = np.sum(adata.X, axis=1)
    tmp_array = adata.X[adata.X != 0] = 1
    sum_peaks = np.sum(tmp_array, axis=1)
    if len(sum_peaks) == 1:
        sum_peaks = sum_peaks.tolist()
        sum_peaks = [item for sublist in sum_peaks for item in sublist]
    adata.obs[key_added] = sum_peaks
    
    # number of peaks in a cell
    np.histogram(adata.obs[key_added])
    
    
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
        plt.xlabel('number of features (log scale)')
        fig = plt.hist(np.log(adata.obs[key_added]), bins, color=color, edgecolor=edgecolor)
    else:
        fig = plt.hist(adata.obs[key_added], bins, color=color, edgecolor=edgecolor)
    
    #fig = plot.get_figure()
    if save!= None:
        plt.savefig(save)
plt.show()
    
    
def binarize(adata, copy=False):
    """convert into a binary matrix"""

    if copy:
        adata2 = adata.copy()
        adata2.X[adata2.X != 0] = 1
        return(adata2)
    else:
        adata.X[adata.X != 0] = 1

        

def cal_var(adata, show=True):
    
    if str(type(adata.X)) != "numpy.ndarray":
        adata.var['n_cells'] = adata.X.sum(axis=0)
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = abs(adata.var['prop_shared_cells']-0.5)
    else:
        adata.var['n_cells'] = adata.X.sum(axis=0).tolist()[0]
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

def select_var_feature(adata, min_score=0.5, nb_features=None, show=True, inplace=True):
    """
    This function compute a score to rank the most shared features across all cells. 
    Then it selects the most variable features according to a minimum variance or select a given number of features.
    
    adata: adata object
    min_score: minimum threshold to retain features
    nb_features: default value is None, if specify it will select a the top most variable features.
    if the nb_features is larger than the total number of feature, it filters based on the min_score argument
    show: default value True, it will plot the distribution of var.
    copy: overwrite the adata object or return another object
    """
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