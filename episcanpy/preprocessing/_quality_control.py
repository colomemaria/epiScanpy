import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

def coverage_cells(adata, bins=50, key_added=None):
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
    plt.hist(adata.obs[key_added], bins)
    plt.show
    

def commoness_features(adata, threshold=None, bw=0.5, key_added=None):
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
    sns.kdeplot(np.array(adata.var[key_added]), bw=bw_param)
    
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
        
