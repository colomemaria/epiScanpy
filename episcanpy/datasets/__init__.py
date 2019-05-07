import anndata as ad
import scanpy.api as sc
from anndata import AnnData
import numpy as np


def normal_blobs(n_variables=11, n_centers=5, cluster_std=1.0, n_observations=640) -> AnnData:
    """Gaussian Blobs, ATAC.
    Parameters
    ----------
    n_variables : `int`, optional (default: 11)
        Dimension of feature space.
    n_centers : `int`, optional (default: 5)
        Number of cluster centers.
    cluster_std : `float`, optional (default: 1.0)
        Standard deviation of clusters.
    n_observations : `int`, optional (default: 640)
        Number of observations. By default, this is the same observation number as in
        ``sc.datasets.krumsiek11()``.
    Returns
    -------
    Annotated data matrix containing a observation annotation 'blobs' that
    indicates cluster identity.
    """
    import sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=n_observations,
                                       n_features=n_variables,
                                       centers=n_centers,
                                       cluster_std=cluster_std,
                                       random_state=0)
    return AnnData(X, obs={'blobs': y.astype(str)})


def bad_atac_sim(n_variables=5000, n_centers=4, cluster_std=2.0, n_observations=3000, dropout=None):
    """
    Generate isotropic Gaussian blobs from sklearn 'make_blobs' for clustering. The data are then binarised
    and we account for dropout. 
    n_variables = nb of features for the simulation
    n_centers = number of predicted clusters
    cluster_std = standard deviation of the clusters
    n_observations = number of cells
    dropouts = amount of dropout (between 0 and 1) default=0.5
    
    return a not annotated adata mastrix
    """
    import sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=n_observations,
                                       n_features=n_variables,
                                       centers=n_centers,
                                       center_box = (0, 1),
                                       cluster_std=cluster_std,
                                       random_state=0)
    adata = AnnData(X, obs={'blobs': y.astype(str)})

    if (dropout==None) or (dropout == 0.5):
        threshold = np.median(adata.X)  
    elif dropout>=1:
        threshold = np.amax(adata.X)
    elif dropout>=0.5:
        threshold = np.amax(adata.X)*dropout
    elif dropout<=0:
        threshold = 0
    else:
        threshold = np.median(adata.X)*dropout
    
    adata.X = (adata.X > threshold).astype(np.int_)

    return (adata)

#############################################
def sim_atac(nb_cells=None, nb_features=None, nb_cell_type=None):
    """Generate a simulated ATACseq count matrix (peaks)
    
    Returns
    -------
    ATACseq binary data matrix.
    """
    mport sklearn.datasets
    X, y = sklearn.datasets.make_blobs(n_samples=n_observations,
                                       n_features=n_variables,
                                       centers=n_centers,
                                       cluster_std=cluster_std,
                                       random_state=0)
    return AnnData(X, obs={'blobs': y.astype(str)})

def sim_met(nb_cells=None, nb_features=None, nb_cell_type=None):
    """Generate a simulated CG methylation count matrix (promoters)
    
    Returns
    -------
    Methylation data matrix.
    """
    [...]
    return(count_matrix)

###############################################
def Luo17_raw():
    """Methylation summaries for 5 cells from the mouse prefrontal cortex [Luo17]_.
    
    Returns
    -------
    Methylation summaries name files.
    """
    [...]
    return(methylation_summaries)

def Luo17promoters_raw():
    """CG Methylation levels at promoters. Mouse prefrontal cortex neurons [Luo17]_.
    
    Returns
    -------
    Annotated data matrix (not filtered).
    """ 
    ## Here I put the raw count matrix. Not the annotated/filtered one 
  [...]
  return(count_matrix)
  
def Luo17promoters():
    """CG Methylation levels at promoters. Mouse prefrontal cortex neurons [Luo17]_.
    
    Returns
    -------
    Annotated data matrix.
    """ 
    return(count_matrix)
 # I will need to add the hematopoitic count matrix from Cusanovich and the peak count matrix for brain data
 # Also Cusanovich et al. 2018
