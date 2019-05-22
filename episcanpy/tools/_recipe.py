import scanpy.api as sc
import anndata as ad

def lazy(adata, pp_pca=True, copy=False):
    '''
    Automatically computes PCA coordinates, loadings and variance decomposition, a neighborhood graph of observations,
    t-distributed stochastic neighborhood embedding (tSNE) Uniform Manifold Approximation and Projection (UMAP)
    
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    pp_pca :  `bool` (default: `True`)
        Computes PCA coordinates before the neighborhood graph
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.
    
    
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    X_pca : `adata.obsm`
        PCA coordinates of data.
    connectivities : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    distances : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    X_tsne : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    X_umap : `adata.obsm`
        UMAP coordinates of data.
    
    
    '''
    if copy:
        adata = adata.copy() 
    else:
        adata

    if pp_pca:
        sc.pp.pca(adata)
    
    sc.pp.neighbors(adata)
    sc.tl.pca(adata)
    sc.tl.tsne(adata)
    sc.tl.umap(adata)
    
    if copy:
        return(adata)
    else:
        None
