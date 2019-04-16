import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

def coverage_cells(adata, bins=50, key_added=None):
    if key_added == None:
        key_added='sum_peaks'
    # make sum peaks
    sum_peaks = []
    matrix = np.matrix(adata.X)
    matrix = matrix.tolist()
    for var in matrix:
        sum_peaks.append(sum(var))
    adata.obs[key_added] = sum_peaks
    
    # number of peaks in a cell
    np.histogram(adata.obs[key_added])
    plt.hist(adata.obs[key_added], bins)
    plt.show
    
def commoness_features(adata, threshold=None, bw=0.5, key_added=None):
    if key_added == None:
        key_added='commonness'
    
    adata_tmp = adata.copy()
    adata_tmp = adata_tmp.transpose()
    common = []
    matrix = np.matrix(adata_tmp.X)
    matrix = matrix.transpose()
    matrix = matrix.tolist()
    for var in matrix:
        common.append(sum(var))
    adata.var[key_added] = common

    if threshold != None:
        plt.axvline(x=threshold, color='r')
    bw_param = bw
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(adata.var[key_added]), bw=bw_param)
    
    
def binarize(adata, copy=False):
    # convert into a binary matrix
    threshold, upper, lower = 1.0, 1.0, 0.0
    admatrix = adata.X
    admatrix = np.where(admatrix>threshold, upper, lower)
    if copy:
        adata2 = adata.copy()
        adata2.X = admatrix
        return(adata2)
    else:
        adata.X = admatrix
        
