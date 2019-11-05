import scanpy as sc
from anndata import AnnData
from typing import Optional, Union  # Special
from typing import Sequence, Collection, Iterable  # ABCs
from typing import Tuple, List  # Classes


## everything here is copied from scanpy

# pca, diffmap, draw_graph, tsne, umap --> I also need it in tools
# I also need louvain and violin and matrixplot and heatmap and rank_gene_groups version
# and leiden 


## ADD DPT

def pca(adata, n_comps=50, zero_center=True, svd_solver='auto', random_state=0,
	return_info=False, use_highly_variable=False, dtype='float32', copy=False,
	chunked=False, chunk_size=None):
    """
    See epi.pp.pca
    """
    if copy:
        return(sc.tl.pca(adata, n_comps, zero_center, svd_solver, random_state,
            return_info, use_highly_variable, dtype, copy, chunked, chunk_size))
    else:
        sc.tl.pca(adata, n_comps, zero_center, svd_solver, random_state,
            return_info, use_highly_variable, dtype, copy, chunked, chunk_size)



def diffmap(adata, n_comps=15, copy=False):
    """Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf18]_.

    Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell
    data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested
    by [Haghverdi16]_ in the implementation of [Wolf18]_.

    The width ("sigma") of the connectivity kernel is implicitly determined by
    the number of neighbors used to compute the single-cell graph in
    :func:`~scanpy.pp.neighbors`. To reproduce the original implementation
    using a Gaussian kernel, use `method=='gauss'` in
    :func:`~scanpy.pp.neighbors`. To use an exponential kernel, use the default
    `method=='umap'`. Differences between these options shouldn't usually be
    dramatic.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_comps : `int`, optional (default: 15)
        The number of dimensions of the representation.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_diffmap** : :class:`numpy.ndarray` (`adata.obsm`)
        Diffusion map representation of data, which is the right eigen basis of
        the transition matrix with eigenvectors as columns.
    **diffmap_evals** : :class:`numpy.ndarray` (`adata.uns`)
        Array of size (number of eigen vectors). Eigenvalues of transition matrix.
    """
    if copy:
        return(sc.tl.diffmap(adata, n_comps, copy))
    else:
        sc.tl.diffmap(adata, n_comps, copy)


def draw_graph(adata, layout= 'fa', init_pos = None, root = None,
    random_state = 0, n_jobs = None, adjacency = None, key_added_ext = None,
    copy = False):

    """\
    Force-directed graph drawing [Islam11]_ [Jacomy14]_ [Chippada18]_.

    An alternative to tSNE that often preserves the topology of the data
    better. This requires to run :func:`~scanpy.pp.neighbors`, first.

    The default layout ('fa', `ForceAtlas2`) [Jacomy14]_ uses the package `fa2
    <https://github.com/bhargavchippada/forceatlas2>`__ [Chippada18]_, which can
    be installed via `pip install fa2`.

    `Force-directed graph drawing
    <https://en.wikipedia.org/wiki/Force-directed_graph_drawing>`__ describes a
    class of long-established algorithms for visualizing graphs. It has been
    suggested for visualizing single-cell data by [Islam11]_. Many other layouts
    as implemented in igraph [Csardi06]_ are available. Similar approaches have
    been used by [Zunder15]_ or [Weinreb17]_.

    Parameters
    ----------
    adata
        Annotated data matrix.
    layout
        'fa' (`ForceAtlas2`) or any valid `igraph layout
        <http://igraph.org/c/doc/igraph-Layout.html>`__. Of particular interest
        are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid Fruchterman Reingold,
        faster than 'fr'), 'kk' (Kamadi Kawai', slower than 'fr'), 'lgl' (Large
        Graph, very fast), 'drl' (Distributed Recursive Layout, pretty fast) and
        'rt' (Reingold Tilford tree layout).
    root
        Root for tree layouts.
    random_state
        For layouts with random initialization like 'fr', change this to use
        different intial states for the optimization. If `None`, no seed is set.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']`.
    key_added_ext
        By default, append `layout`.
    proceed
        Continue computation, starting off with 'X_draw_graph_`layout`'.
    init_pos
        `'paga'`/`True`, `None`/`False`, or any valid 2d-`.obsm` key.
        Use precomputed coordinates for initialization.
        If `False`/`None` (the default), initialize randomly.
    copy
        Return a copy instead of writing to adata.
    **kwds
        Parameters of chosen igraph layout. See e.g. `fruchterman-reingold`_
        [Fruchterman91]_. One of the most important ones is `maxiter`.

        .. _fruchterman-reingold: http://igraph.org/python/doc/igraph.Graph-class.html#layout_fruchterman_reingold

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following field.

    **X_draw_graph_layout** : `adata.obsm`
        Coordinates of graph layout. E.g. for layout='fa' (the default), the field is called 'X_draw_graph_fa'
    """
    if copy:
        return(sc.tl.draw_graph(adata, layout, init_pos, root, random_state, n_jobs, 
            adjacency, key_added_ext, copy))
    else:
        sc.tl.draw_graph(adata, layout, init_pos, root, random_state, n_jobs,
            adjacency, key_added_ext, copy)

def tsne(adata,
    n_pcs=None,
    use_rep=None,
    perplexity=30,
    early_exaggeration=12,
    learning_rate=1000,
    random_state=0,
    use_fast_tsne=True,
    n_jobs=None,
    copy=False):
    """\
    t-SNE [Maaten08]_ [Amir13]_ [Pedregosa11]_.

    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08]_ has been
    proposed for visualizating single-cell data by [Amir13]_. Here, by default,
    we use the implementation of *scikit-learn* [Pedregosa11]_. You can achieve
    a huge speedup and better convergence if you install `Multicore-tSNE
    <https://github.com/DmitryUlyanov/Multicore-TSNE>`__ by [Ulyanov16]_, which
    will be automatically detected by Scanpy.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    {doc_n_pcs}
    {use_rep}
    perplexity : `float`, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.
    early_exaggeration : `float`, optional (default: 12.0)
        Controls how tight natural clusters in the original space are in the
        embedded space and how much space will be between them. For larger
        values, the space between natural clusters will be larger in the
        embedded space. Again, the choice of this parameter is not very
        critical. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate : `float`, optional (default: 1000)
        Note that the R-package "Rtsne" uses a default of 200.
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.
    random_state : `int` or `None`, optional (default: 0)
        Change this to use different intial states for the optimization. If `None`,
        the initial state is not reproducible.
    use_fast_tsne : `bool`, optional (default: `True`)
        Use the MulticoreTSNE package by D. Ulyanov if it is installed.
    n_jobs : `int` or `None` (default: `sc.settings.n_jobs`)
        Number of jobs.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_tsne** : `np.ndarray` (`adata.obs`, dtype `float`)
        tSNE coordinates of data.
    """
    if copy:
        return(sc.tl.tsne(adata, n_pcs, use_rep, perplexity, early_exaggeration, 
            learning_rate, random_state, use_fast_tsne, n_jobs, copy))
    else:
        sc.tl.tsne(adata, n_pcs, use_rep, perplexity, early_exaggeration,
            learning_rate, random_state, use_fast_tsne, n_jobs, copy)


def umap(adata,
    min_dist=0.5,
    spread=1.0,
    n_components=2,
    maxiter=None,
    alpha=1.0,
    gamma=1.0,
    negative_sample_rate=5,
    init_pos='spectral',
    random_state=0,
    a=None,
    b=None,
    copy=False,
    #method='gauss',
    ):

    """Embed the neighborhood graph using UMAP [McInnes18]_.

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.  We use the
    implementation of `umap-learn <https://github.com/lmcinnes/umap>`__
    [McInnes18]_. For a few comparisons of UMAP with tSNE, see this `preprint
    <https://doi.org/10.1101/298430>`__.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    min_dist : `float`, optional (default: 0.5)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points on
        the manifold are drawn closer together, while larger values will result
        on a more even dispersal of points. The value should be set relative to
        the ``spread`` value, which determines the scale at which embedded
        points will be spread out. The default of in the `umap-learn` package is
        0.1.
    spread : `float` (optional, default 1.0)
        The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are.
    n_components : `int`, optional (default: 2)
        The number of dimensions of the embedding.
    maxiter : `int`, optional (default: `None`)
        The number of iterations (epochs) of the optimization. Called `n_epochs`
        in the original UMAP.
    alpha : `float`, optional (default: 1.0)
        The initial learning rate for the embedding optimization.
    gamma : `float` (optional, default 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate : `int` (optional, default 5)
        The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding.
    init_pos : `string` or `np.array`, optional (default: 'spectral')
        How to initialize the low dimensional embedding. Called `init` in the
        original UMAP.
        Options are:

        * Any key for `adata.obsm`.
        * 'paga': positions from :func:`~scanpy.pl.paga`.
        * 'spectral': use a spectral embedding of the graph.
        * 'random': assign initial embedding positions at random.
        * A numpy array of initial embedding positions.
    random_state : `int`, `RandomState` or `None`, optional (default: 0)
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a : `float` (optional, default `None`)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    b : `float` (optional, default `None`)
        More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.
    method : {`'umap'`, `'rapids'`}  (default: `'gauss'`)
        Use the original 'umap' implementation, or 'rapids' (experimental, GPU only)

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_umap** : `adata.obsm` field
        UMAP coordinates of data.
    """
    if copy:
        return(sc.tl.umap(adata, min_dist, spread, n_components, maxiter, alpha, gamma,
            negative_sample_rate, init_pos, random_state, a, b,copy))
    else:
        sc.tl.umap(adata, min_dist, spread, n_components, maxiter, alpha, gamma,
            negative_sample_rate, init_pos, random_state, a, b,copy)

def louvain(adata, resolution=None, random_state=0, restrict_to=None, 
	key_added='louvain', adjacency=None, flavor='vtraag', directed=True,
	use_weights=False, partition_type=None, partition_kwargs=None, copy=False):
    """Cluster cells into subgroups [Blondel08]_ [Levine15]_ [Traag17]_.

    Cluster cells using the Louvain algorithm [Blondel08]_ in the implementation
    of [Traag17]_. The Louvain algorithm has been proposed for single-cell
    analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or :func:`~scanpy.external.pp.bbknn` first,
    or explicitly passing a ``adjacency`` matrix.

    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        For the default flavor (``'vtraag'``), you can provide a resolution
        (higher resolution means finding more and smaller clusters),
        which defaults to 1.0. See “Time as a resolution parameter” in [Lambiotte09]_.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain ``(obs_key, list_of_categories)``.
    key_added
        Key under which to add the cluster labels. (default: ``'louvain'``)
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        ``adata.uns['neighbors']['connectivities']``.
    flavor : {``'vtraag'``, ``'igraph'``, ``'rapids'``}
        Choose between to packages for computing the clustering.
        ``'vtraag'`` is much more powerful, and the default.
    directed
        Interpret the ``adjacency`` matrix as directed graph?
    use_weights
        Use weights from knn graph.
    partition_type
        Type of partition to use.
        Only a valid argument if ``flavor`` is ``'vtraag'``.
    partition_kwargs
        Key word arguments to pass to partitioning,
        if ``vtraag`` method is being used.
    copy
        Copy adata or modify it inplace.

    Returns
    -------
    :obj:`None`
        By default (``copy=False``), updates ``adata`` with the following fields:

        ``adata.obs['louvain']`` (:class:`pandas.Series`, dtype ``category``)
            Array of dim (number of samples) that stores the subgroup id
            (``'0'``, ``'1'``, ...) for each cell.

    :class:`~anndata.AnnData`
        When ``copy=True`` is set, a copy of ``adata`` with those fields is returned.
    """

    if copy:
        return(sc.tl.louvain(adata, resolution, random_state, restrict_to,
            key_added, adjacency, flavor, directed, use_weights,
            partition_type, partition_kwargs, copy))
    else:
        sc.tl.louvain(adata, resolution, random_state, restrict_to,
            key_added, adjacency, flavor, directed, use_weights,
            partition_type, partition_kwargs, copy)

def leiden(adata,
    resolution= 1,
    *,
    restrict_to = None,
    random_state = 0,
    key_added = 'leiden',
    adjacency = None,
    directed = True,
    use_weights = True,
    n_iterations = -1,
    partition_type = None,
    copy = False):

    """Cluster cells into subgroups [Traag18]_.

    Cluster cells using the Leiden algorithm [Traag18]_, an improved version of the
    Louvain algorithm [Blondel08]_. The Louvain algorithm has been proposed for single-cell
    analysis by [Levine15]_.

    This requires having ran :func:`~scanpy.pp.neighbors` or :func:`~scanpy.external.pp.bbknn` first.

    Parameters
    ----------
    adata
        The annotated data matrix.
    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters. Set to `None` if overriding `partition_type`
        to one that doesn’t accept a `resolution_parameter`.
    random_state
        Change the initialization of the optimization.
    restrict_to
        Restrict the clustering to the categories within the key for sample
        annotation, tuple needs to contain `(obs_key, list_of_categories)`.
    key_added
        `adata.obs` key under which to add the cluster labels. (default: `'leiden'`)
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        `adata.uns['neighbors']['connectivities']`.
    directed
        Whether to treat the graph as directed or undirected.
    use_weights
        If `True`, edge weights from the graph are used in the computation
        (placing more emphasis on stronger edges).
    n_iterations
        How many iterations of the Leiden clustering algorithm to perform.
        Positive values above 2 define the total number of iterations to perform,
        -1 has the algorithm run until it reaches its optimal clustering.
    partition_type
        Type of partition to use. Defaults to :class:`~leidenalg.RBConfigurationVertexPartition`.
        For the available options, consult the documentation for :func:`~leidenalg.find_partition`.
    copy
        Whether to copy `adata` or modify it inplace.
    **partition_kwargs
        Any further arguments to pass to `~leidenalg.find_partition`
        (which in turn passes arguments to the `partition_type`).

    Returns
    -------
    `adata.obs[key_added]`
        Array of dim (number of samples) that stores the subgroup id (`'0'`, `'1'`, ...) for each cell.
    `adata.uns['leiden']['params']`
        A dict with the values for the parameters `resolution`, `random_state`, and `n_iterations`.
    """

    if copy:
        return(sc.tl.leiden(adata, resolution, #restrict_to, random_state, key_added, 
            #adjacency, directed, use_weights, n_iterations, partition_type, copy
            ))
    else:
        sc.tl.leiden(adata, resolution, #restrict_to, random_state, key_added, 
            #adjacency, directed, use_weights, n_iterations, partition_type, copy
            )


def dendogram(
    adata: AnnData,
    groupby: str,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    var_names: Optional[Sequence[str]] = None,
    use_raw: Optional[bool] = None,
    cor_method: str = 'pearson',
    linkage_method: str = 'complete',
    key_added: Optional[str] = None,
) -> None:
    """\
    Computes a hierarchical clustering for the given `groupby` categories.

    By default, the PCA representation is used unless `.X` has less than 50 variables.

    Alternatively, a list of `var_names` (e.g. genes) can be given.

    Average values of either `var_names` or components are used to compute a correlation matrix.

    The hierarchical clustering can be visualized using `sc.pl.dendrogram` or multiple other
    visualizations that can include a dendrogram: `matrixplot`, `heatmap`, `dotplot` and `stacked_violin`

    .. note::
        The computation of the hierarchical clustering is based on predefined groups and not
        per cell. The correlation matrix is computed using by default pearson but other methods
        are available.

    Parameters
    ----------
    adata
        Annotated data matrix
    {n_pcs}
    {use_rep}
    var_names
        List of var_names to use for computing the hierarchical clustering.
        If `var_names` is given, then `use_rep` and `n_pcs` is ignored.
    use_raw
        Only when `var_names` is not None.
        Use `raw` attribute of `adata` if present.
    cor_method
        correlation method to use.
        Options are 'pearson', 'kendall', and 'spearman'
    linkage_method
        linkage method to use. See :func:`scipy.cluster.hierarchy.linkage`
        for more information.
    key_added
        By default, the dendrogram information is added to
        `.uns[f'dendrogram_{{groupby}}']`.
        Notice that the `groupby` information is added to the dendrogram.

    Returns
    -------
    `adata.uns['dendrogram']` (or instead of 'dendrogram' the value selected
    for `key_added`) is updated with the dendrogram information

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, groupby='bulk_labels')
    >>> sc.pl.dendrogram(adata)
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)
    """

    sc.tl.dendrogram(adata=adata,
    groupby=groupby,
    n_pcs=n_pcs,
    use_rep=use_rep,
    var_names=var_names,
    use_raw=use_raw,
    cor_method=cor_method,
    linkage_method=linkage_method,
    key_added=key_added)


