import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import os
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import homogeneity_score
import seaborn as sns


def getNClusters(adata,n_cluster,range_min=0,range_max=3,max_steps=20, method='louvain', key_added=None):
    """
    Function will test different settings of louvain to obtain the target number of clusters.
    credit to Lucas Pinello lab. See: https://github.com/pinellolab/scATAC-benchmarking

    It can get cluster for both louvain and leiden.
    You can specify the obs variable name as key_added. 
    """
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        print('step ' + str(this_step))
        this_resolution = this_min + ((this_max-this_min)/2)
        
        if (method == 'louvain') and (key_added==None):
            sc.tl.louvain(adata, resolution=this_resolution)
        elif method == 'louvain':
            sc.tl.louvain(adata, resolution=this_resolution, key_added=key_added)
        elif( method == 'leiden') and (key_added==None):
            sc.tl.leiden(adata,resolution=this_resolution)
        else:
            sc.tl.leiden(adata,resolution=this_resolution, key_added=key_added)
    
        
        if key_added==None:
            this_clusters = adata.obs[method].nunique()
        else:
            this_clusters = adata.obs[key_added].nunique()
        
        print('got ' + str(this_clusters) + ' at resolution ' + str(this_resolution))
        
        if this_clusters > n_cluster:
            this_max = this_resolution
        elif this_clusters < n_cluster:
            this_min = this_resolution
        else:
            return(this_resolution, adata)
        this_step += 1
    
    print('Cannot find the number of clusters')
    print('Clustering solution from last iteration is used:' + str(this_clusters) + ' at resolution ' + str(this_resolution))
   

def kmeans(adata, num_clusters):
    """
    Compute kmeans clustering using X_pca fits.
    random_state = 2019
    """
    kmeans = KMeans(n_clusters=num_clusters, random_state=2019).fit(adata.obsm['X_pca']) 
    adata.obs['kmeans'] = pd.Series(kmeans.labels_,index=adata.obs.index).astype('category')

def hc(adata, num_clusters):
    """
    Compute hierarchical clustering using X_pca fits.
    random_state = 2019
    """
    hc = AgglomerativeClustering(n_clusters=num_clusters).fit(adata.obsm['X_pca'])
    adata.obs['hc'] = pd.Series(hc.labels_,index=adata.obs.index).astype('category')

#### Metrics for clustering

def ARI(adata, label_1, label_2):
    """
    Compute Adjusted Rand Index.
    """

    return(adjusted_rand_score(adata.obs[label_1], adata.obs[label_2]))


def AMI(adata, label_1, label_2):
    """
    Compute adjusted Mutual Info.
    """

    return(adjusted_mutual_info_score(adata.obs[label_1], adata.obs[label_2]))

def homogeneity(adata, label_1, label_2):
    """
    Compute homogeneity score.
    """

    return(homogeneity_score(adata.obs[label_1], adata.obs[label_2]))

