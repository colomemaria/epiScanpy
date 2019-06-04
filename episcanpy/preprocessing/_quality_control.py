import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

def commonness_features(adata, threshold=None, bw=0.5, key_added=None, xlabel=None, title=None, save=None):
    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.
    """
    if key_added == None:
        key_added='commonness'
    
    common = np.sum(adata.X, axis=0).tolist()
    if len(common) == 1:
        common = [item for sublist in common for item in sublist]
    adata.var[key_added] = common

    if threshold != None:
        plt.axvline(x=threshold, color='r')
    bw_param = bw
    sns.set_style('whitegrid')
    #sns.kdeplot(np.array(adata.var[key_added]), bw=bw_param)
    
    plot = sns.distplot(common, hist=False, kde=True,
                    bins=int(80), color = 'darkblue', 
                    hist_kws={'edgecolor':'black'},
                    kde_kws={'linewidth': 1})
    if xlabel ==None:
        plt.xlabel('cells sharing a feature')
    else:
        plt.xlabel(xlabel)
        
    
    if title !=None:
        plt.title(title)
    
    fig = plot.get_figure()
    if save!= None:
        fig.savefig(save)
    plt.show()
    
def coverage_cells(adata, bins=50, key_added=None, xlabel=None, ylabel=None, title=None, save=None):
    """
    Histogram of the number of features with an open peak (in the case of ATAC-seq data)
    for every cells.
    """
    if key_added == None:
        key_added='sum_peaks'
    # make sum peaks
    
    sum_peaks = np.sum(adata.X, axis=1)
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
        
    fig = plt.hist(adata.obs[key_added], bins)
    
    #fig = plot.get_figure()
    if save!= None:
        plt.savefig(save)
    plt.show()
    
    
def binarize(adata, copy=False):
    """convert into a binary matrix"""
    threshold, upper, lower = 1.0, 1.0, 0.0
    admatrix = adata.X
    admatrix = np.where(admatrix>threshold, upper, lower)
    if copy:
        adata2 = adata.copy()
        adata2.X = admatrix
        return(adata2)
    else:
        adata.X = admatrix
        
