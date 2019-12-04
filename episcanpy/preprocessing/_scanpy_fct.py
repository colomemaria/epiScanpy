## all functions here are a direct copy from Scanpy version ### (October 18)
import scanpy as sc

N_PCS = 50

def filter_cells(
    adata,
    min_counts = None,
    min_features = None,
    max_counts = None,
    max_features = None,
    inplace = True,
    copy = False):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_features` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters ``min_counts``, ``min_features``,
    ``max_counts``, ``max_features`` per call.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    min_counts
        Minimum number of counts required for a cell to pass filtering.
    min_features
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts
        Maximum number of counts required for a cell to pass filtering.
    max_features
        Maximum number of genes expressed required for a cell to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on ``inplace``, returns the following arrays or directly subsets
    and annotates the data matrix:

    cells_subset : numpy.ndarray
        Boolean index mask that does filtering. ``True`` means that the
        cell is kept. ``False`` means the cell is removed.
    number_per_cell : numpy.ndarray
        Depending on what was tresholded (``counts`` or ``genes``), the array stores
        ``n_counts`` or ``n_cells`` per gene.

    Examples
    --------
    >>> adata = sc.datasets.krumsiek11()
    >>> adata.n_obs
    640
    >>> adata.var_names
    ['Gata2' 'Gata1' 'Fog1' 'EKLF' 'Fli1' 'SCL' 'Cebpa'
     'Pu.1' 'cJun' 'EgrNab' 'Gfi1']
    >>> # add some true zeros
    >>> adata.X[adata.X < 0.3] = 0
    >>> # simply compute the number of genes per cell
    >>> sc.pp.filter_cells(adata, min_features=0)
    >>> adata.n_obs
    640
    >>> adata.obs['nb_features'].min()
    1
    >>> # filter manually
    >>> adata_copy = adata[adata.obs['nb_features'] >= 3]
    >>> adata_copy.obs['nb_features'].min()
    >>> adata.n_obs
    554
    >>> adata.obs['nb_features'].min()
    3
    >>> # actually do some filtering
    >>> sc.pp.filter_cells(adata, min_features=3)
    >>> adata.n_obs
    554
    >>> adata.obs['nb_features'].min()
    3
    """
    if copy:
        adata_copy = sc.pp.filter_cells(adata, min_counts, min_genes=min_features, max_counts=max_counts,
            max_genes=max_features, inplace=inplace, copy=copy)
        adata_copy.obs['nb_features'] = adata_copy.obs['n_genes'] 
        del adata_copy.obs['n_genes']
        return(adata_copy)
    else:
        sc.pp.filter_cells(adata, min_counts, min_genes=min_features, max_counts=max_counts, 
            max_genes=max_features, inplace=inplace, copy=copy)
        adata.obs['nb_features'] = adata.obs['n_genes'] 
        del adata.obs['n_genes']



def filter_features(data,
    min_counts = None,
    min_cells = None,
    max_counts = None,
    max_cells = None,
    inplace = True,
    copy = False):
    """Filter features based on number of cells or counts.

    Keep features that have at least ``min_counts`` counts or are expressed in at
    least ``min_cells`` cells or have at most ``max_counts`` counts or are expressed
    in at most ``max_cells`` cells.

    Only provide one of the optional parameters ``min_counts``, ``min_cells``,
    ``max_counts``, ``max_cells`` per call.

    Parameters
    ----------
    data
        An annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts
        Minimum number of counts required for a gene to pass filtering.
    min_cells
        Minimum number of cells expressed required for a gene to pass filtering.
    max_counts
        Maximum number of counts required for a gene to pass filtering.
    max_cells
        Maximum number of cells expressed required for a gene to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on `inplace`, returns the following arrays or directly subsets
    and annotates the data matrix

    gene_subset : numpy.ndarray
        Boolean index mask that does filtering. `True` means that the
        gene is kept. `False` means the gene is removed.
    number_per_gene : numpy.ndarray
        Depending on what was tresholded (`counts` or `cells`), the array stores
        `n_counts` or `n_cells` per gene.
    """
    if copy:
        return(sc.pp.filter_genes(data, min_counts, min_cells, max_counts,
            max_cells, inplace, copy))
    else:
        sc.pp.filter_genes(data, min_counts, min_cells, max_counts, max_cells, inplace, copy)
    


def pca(adata,
    n_comps = N_PCS,
    zero_center = True,
    svd_solver = 'auto',
    random_state = 0,
    return_info = False,
    use_highly_variable = False,
    dtype = 'float32',
    copy = False,
    chunked = False,
    chunk_size = None):
    """Principal component analysis [Pedregosa11]_.

    Computes PCA coordinates, loadings and variance decomposition. Uses the
    implementation of *scikit-learn* [Pedregosa11]_.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute.
    zero_center
        If `True`, compute standard PCA from covariance matrix.
        If ``False``, omit zero-centering variables
        (uses :class:`~sklearn.decomposition.TruncatedSVD`),
        which allows to handle sparse input efficiently.
        Passing ``None`` decides automatically based on sparseness of the data.
    svd_solver
        SVD solver to use:

        ``'arpack'``
          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)

        ``'randomized'``
          for the randomized algorithm due to Halko (2009).

        ``'auto'`` (the default)
          chooses automatically depending on the size of the problem.

    random_state
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        ``.var['highly_variable']``.
        By default uses them if they have been determined beforehand.
    dtype
        Numpy data type string to which to convert the result.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    chunked
        If ``True``, perform an incremental PCA on segments of ``chunk_size``.
        The incremental PCA automatically zero centers and ignores settings of
        ``random_seed`` and ``svd_solver``. If ``False``, perform a full PCA.
    chunk_size
        Number of observations to include in each chunk.
        Required if ``chunked=True`` was passed.

    Returns
    -------
    X_pca : :class:`scipy.sparse.spmatrix` or :class:`numpy.ndarray`
        If `data` is array-like and ``return_info=False`` was passed,
        this function only returns `X_pca`…
    adata : anndata.AnnData
        …otherwise if ``copy=True`` it returns or else adds fields to ``adata``:

        ``.obsm['X_pca']``
             PCA representation of data.

        ``.varm['PCs']``
             The principal components containing the loadings.

        ``.uns['pca']['variance_ratio']``)
             Ratio of explained variance.

        ``.uns['pca']['variance']``
             Explained variance, equivalent to the eigenvalues of the covariance matrix.
    """
    # chunked calculation is not randomized, anyways

    if copy:
        return(sc.pp.pca(adata, n_comps, zero_center, svd_solver, random_state, return_info,
            use_highly_variable, dtype, copy, chunked, chunk_size))
    else:
        sc.pp.pca(adata, n_comps, zero_center, svd_solver, random_state, return_info, 
            use_highly_variable, dtype, copy, chunked, chunk_size)


def normalize_per_cell(
    adata,
    counts_per_cell_after=None,
    counts_per_cell=None,
    key_n_counts=None,
    copy=False,
    layers=[],
    use_rep=None,
    min_counts=1):
    """Normalize total counts per cell.

    .. warning::
        .. deprecated:: 1.3.7
            Use :func:`~scanpy.pp.normalize_total` instead.
            The new function is equivalent to the present
            function, except that

            * the new function doesn't filter cells based on `min_counts`,
              use :func:`~scanpy.pp.filter_cells` if filtering is needed.
            * some arguments were renamed
            * `copy` is replaced by `inplace`

    Normalize each cell by total counts over all genes, so that every cell has
    the same total count after normalization.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    counts_per_cell_after : `float` or `None`, optional (default: `None`)
        If `None`, after normalization, each cell has a total count equal
        to the median of the *counts_per_cell* before normalization.
    counts_per_cell : `np.array`, optional (default: `None`)
        Precomputed counts per cell.
    key_n_counts : `str`, optional (default: `'n_counts'`)
        Name of the field in `adata.obs` where the total counts per cell are
        stored.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.
    min_counts : `int`, optional (default: 1)
        Cells with counts less than `min_counts` are filtered out during
        normalization.

    Returns
    -------
    Returns or updates `adata` with normalized version of the original
    `adata.X`, depending on `copy`.

    Examples
    --------
    >>> adata = AnnData(
    >>>     data=np.array([[1, 0], [3, 0], [5, 6]]))
    >>> print(adata.X.sum(axis=1))
    [  1.   3.  11.]
    >>> sc.pp.normalize_per_cell(adata)
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts
    0       1.0
    1       3.0
    2      11.0
    [ 3.  3.  3.]
    >>> sc.pp.normalize_per_cell(adata, counts_per_cell_after=1,
    >>>                          key_n_counts='n_counts2')
    >>> print(adata.obs)
    >>> print(adata.X.sum(axis=1))
       n_counts  n_counts2
    0       1.0        3.0
    1       3.0        3.0
    2      11.0        3.0
    [ 1.  1.  1.]
    """

    if copy:
        return(sc.pp.normalize_per_cell(adata, counts_per_cell_after, counts_per_cell,
            key_n_counts,copy, layers, use_rep, min_counts))
    else:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after, counts_per_cell,
            key_n_counts,copy, layers, use_rep, min_counts)

def normalize_total(adata,
    target_sum = None,
    exclude_highly_expressed = False,
    max_fraction = 0.05,
    key_added = None,
    layers = None,
    layer_norm = None,
    inplace = True,
):
    """\
    Normalize counts per cell.

    If choosing ``target_sum=1e6``, this is CPM normalization.

    If ``exclude_highly_expressed=True``, very highly expressed genes are excluded
    from the computation of the normalization factor (size factor) for each
    cell. This is meaningful as these can strongly influence the resulting
    normalized values for all other genes [Weinreb17]_.

    Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
    [Zheng17]_ or SPRING [Weinreb17]_.

    Params
    ------
    adata
        The annotated data matrix of shape ``n_obs`` × ``n_vars``. Rows correspond
        to cells and columns to features.
    target_sum
        If ``None``, after normalization, each observation (cell) has a total count
        equal to the median of total counts for observations (cells)
        before normalization.
    exclude_highly_expressed
        Exclude (very) highly expressed genes for the computation of the
        normalization factor (size factor) for each cell. A gene is considered
        highly expressed, if it has more than ``max_fraction`` of the total counts
        in at least one cell. The not-excluded genes will sum up to
        ``target_sum``.
    max_fraction
        If ``exclude_highly_expressed=True``, consider cells as highly expressed
        that have more counts than ``max_fraction`` of the original total counts
        in at least one cell.
    key_added
        Name of the field in ``adata.obs`` where the normalization factor is
        stored.
    layers
        List of layers to normalize. Set to ``'all'`` to normalize all layers.
    layer_norm
        Specifies how to normalize layers:

        * If `None`, after normalization, for each layer in *layers* each cell
          has a total count equal to the median of the *counts_per_cell* before
          normalization of the layer.
        * If `'after'`, for each layer in *layers* each cell has
          a total count equal to `target_sum`.
        * If `'X'`, for each layer in *layers* each cell has a total count
          equal to the median of total counts for observations (cells) of
          `adata.X` before normalization.

    inplace
        Whether to update ``adata`` or return dictionary with normalized copies of
        ``adata.X`` and ``adata.layers``.

    Returns
    -------
    Returns dictionary with normalized copies of `adata.X` and `adata.layers`
    or updates `adata` with normalized version of the original
    `adata.X` and `adata.layers`, depending on `inplace`.

    Example
    --------
    >>> from anndata import AnnData
    >>> import scanpy as sc
    >>> sc.settings.verbosity = 2
    >>> np.set_printoptions(precision=2)
    >>> adata = AnnData(np.array([[3, 3, 3, 6, 6], [1, 1, 1, 2, 2], [1, 22, 1, 2, 2]]))
    >>> adata.X
    array([[ 3.,  3.,  3.,  6.,  6.],
           [ 1.,  1.,  1.,  2.,  2.],
           [ 1., 22.,  1.,  2.,  2.]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, inplace=False)['X']
    >>> X_norm
    array([[0.14, 0.14, 0.14, 0.29, 0.29],
           [0.14, 0.14, 0.14, 0.29, 0.29],
           [0.04, 0.79, 0.04, 0.07, 0.07]], dtype=float32)
    >>> X_norm = sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=True, max_fraction=0.2, inplace=False)['X']
    The following highly-expressed genes are not considered during normalization factor computation:
    ['1', '3', '4']
    >>> X_norm
    array([[ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5,  0.5,  0.5,  1. ,  1. ],
           [ 0.5, 11. ,  0.5,  1. ,  1. ]], dtype=float32)
    """

    sc.pp.normalize_total(adata, target_sum, exclude_highly_expressed, 
        max_fraction, key_added, layers, layer_norm, inplace)
    

def regress_out(adata, keys, n_jobs=None, copy=False):
    """Regress out unwanted sources of variation.

    Uses simple linear regression. This is inspired by Seurat's `regressOut`
    function in R [Satija15].

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    keys : `str` or list of `str`
        Keys for observation annotation on which to regress on.
    n_jobs : `int` or `None`, optional. If None is given, then the n_jobs seting is used (default: `None`)
        Number of jobs for parallel computation.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates `adata` with the corrected data matrix.
    """
    if copy:
        return(sc.pp.regress_out(adata, keys, n_jobs, copy))
    else:
        sc.pp.regress_out(adata, keys, n_jobs, copy)
    

def subsample(data, fraction=None, n_obs=None, random_state=0, copy=False):
    """Subsample to a fraction of the number of observations.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`, `np.ndarray`, `sp.sparse`
        The (annotated) data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    fraction : `float` in [0, 1] or `None`, optional (default: `None`)
        Subsample to this `fraction` of the number of observations.
    n_obs : `int` or `None`, optional (default: `None`)
        Subsample to this number of observations.
    random_state : `int` or `None`, optional (default: 0)
        Random seed to change subsampling.
    copy : `bool`, optional (default: `False`)
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Returns `X[obs_indices], obs_indices` if data is array-like, otherwise
    subsamples the passed :class:`~anndata.AnnData` (`copy == False`) or
    returns a subsampled copy of it (`copy == True`).
    """
    if copy:
        return(sc.pp.subsample(data, fraction, n_obs, random_state, copy))
    else:
        sc.pp.subsample(data, fraction, n_obs, random_state, copy)
    
def downsample_counts(adata, counts_per_cell = None, total_counts = None,
    random_state = 0, replace = False, copy = False):
    """Downsample counts from count matrix.

    If `counts_per_cell` is specified, each cell will downsampled. If
    `total_counts` is specified, expression matrix will be downsampled to
    contain at most `total_counts`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    counts_per_cell
        Target total counts per cell. If a cell has more than 'counts_per_cell',
        it will be downsampled to this number. Resulting counts can be specified
        on a per cell basis by passing an array.Should be an integer or integer
        ndarray with same length as number of obs.
    total_counts
        Target total counts. If the count matrix has more than `total_counts`
        it will be downsampled to have this number.
    random_state
        Random seed for subsampling.
    replace
        Whether to sample the counts with replacement.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned.

    Returns
    -------
    Depending on `copy` returns or updates an `adata` with downsampled `.X`.
    """

    if copy:
        return(sc.pp.downsample_counts(adata, counts_per_cell, total_counts,
            random_state, replace, copy))
    else:
        sc.pp.downsample_counts(adata, counts_per_cell, total_counts,
            random_state, replace, copy)
    
    



#def zscore_deprecated(X: np.ndarray) -> np.ndarray:
#    """Z-score standardize each variable/gene in X.
#
#    Use `scale` instead.

#    Reference: Weinreb et al. (2017).

#    Parameters
#    ----------
#    X
#        Data matrix. Rows correspond to cells and columns to genes.

#    Returns
#    -------
#    Z-score standardized version of the data matrix.
#    """
#    means = np.tile(np.mean(X, axis=0)[None, :], (X.shape[0], 1))
#    stds = np.tile(np.std(X, axis=0)[None, :], (X.shape[0], 1))
#    return (X - means) / (stds + .0001)


