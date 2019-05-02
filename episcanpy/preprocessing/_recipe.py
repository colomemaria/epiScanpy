import scanpy.api as sc
import anndata as ad

def lazy(adata, pp_pca=True, nb_pcs=50, n_neighbors=15, perplexity=30, 
         method='umap', metric='euclidean', min_dist=0.5, spread=1.0, 
         n_components=2, copy=False):
    '''
    Automatically computes PCA coordinates, loadings and variance decomposition, a neighborhood graph of observations,
    t-distributed stochastic neighborhood embedding (tSNE) Uniform Manifold Approximation and Projection (UMAP)
    
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    pp_pca :  `bool` (default: `True`)
        Computes PCA coordinates before the neighborhood graph
    nb_pcs : `int` (default: 50)
        Number of principal component computed for PCA (and therefore neighbors, tsne and umap)
    n_neighbors : `int` (default: 15)
        The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation.
        Larger values result in more global views of the manifold, while smaller values result in more local data being
        preserved. In general values should be in the range 2 to 100.
    method : `str` (default: ‘umap’)
        Use ‘umap’ or ‘gauss’, kernel for computing connectivities. Gives very similar results.
    metric : `str` (default: ‘euclidean’)
        A known metric’s name or a callable that returns a distance.
    perplexity : `int` (default: 30)
        The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms.
        Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. 
    min_dist : `float` (default: 0.5)
        The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped
        embedding where nearby points on the manifold are drawn closer together, while larger values will result on a
        more even dispersal of points.
    spread : `float` (default: 1.0)
        The effective scale of embedded points. In combination with `min_dist` this determines how clustered/clumped
        the embedded points are.
    n_components : `int` (default: 2)
        The number of dimensions of the `UMAP` embedding.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.
    
    method='umap',
    metric='euclidean',
    min_dist=0.5,
    spread=1.0, 
    n_components=2
    
    
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
        sc.pp.pca(adata, n_comps=nb_pcs)
    
    sc.pp.neighbors(adata,  n_neighbors=n_neighbors, n_pcs=nb_pcs, method=method, metric=metric)
    sc.tl.pca(adata, n_comps=nb_pcs)
    sc.tl.tsne(adata, n_pcs=nb_pcs, perplexity=perplexity)
    sc.tl.umap(adata, min_dist, spread, n_components)
    
    if copy:
        return(adata)
    else:
        None
  
