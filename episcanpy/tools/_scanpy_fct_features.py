import scanpy as sc

## everything here is copied from scanpy and adapted

# pca, diffmap, draw_graph, tsne, umap --> I also need it in tools
# heatmap, violin and matrixplot and heatmap and rank_gene_groups version

def stacked_violinstacked_violin(
    adata: AnnData,
    var_names,
    groupby=None,
    log=False,
    use_raw=None,
    num_categories=7,
    figsize=None,
    dendrogram=False,
    feature_symbols=None,
    var_group_positions=None,
    var_group_labels=None,
    standard_scale: Optional[str] = None,
    var_group_rotation=None,
    layer=None,
    stripplot: bool = False,
    jitter: Union[float, bool] = False,
    size: int = 1,
    scale: str = 'width',
    order: Optional[Sequence[str]] = None,
    swap_axes: bool = False,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    row_palette: str = 'muted',
    **kwds,
):
	"""\
    Stacked violin plots.

    Makes a compact image composed of individual violin plots
    (from :func:`~seaborn.violinplot`) stacked on top of each other.
    Useful to visualize gene expression per cluster.

    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.

    Parameters
    ----------
    {common_plot_args}
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    order
        Order in which to show the categories.
    scale: {{`'area'`, `'count'`, `'width'`}}
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    row_palette
        The row palette determines the colors to use for the stacked violins.
        The value should be a valid seaborn or matplotlib palette name
        (see :func:`~seaborn.color_palette`).
        Alternatively, a single color name or hex value can be passed,
        e.g. `'red'` or `'#cc33ff'`.
    standard_scale: {{`'var'`, `'obs'`}}
        Whether or not to standardize a dimension between 0 and 1,
        meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby` categories.
         By setting `swap_axes` then x are the `groupby` categories and y the `var_names`. When swapping
         axes var_group_positions are no longer used
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.violinplot`.

    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`

    Examples
    -------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)

    Using var_names as dict:

    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)

    See also
    --------
    rank_genes_groups_stacked_violin: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """

    sc.pl.stacked_violin(stacked_violin(adata, 
    	var_names=var_names,
    	groupby=groupby,
    	log=log,
    	use_raw=use_raw,
    	num_categories=num_categories,
    	figsize=figsize,
    	dendrogram=dendrogram,
    	gene_symbols=feature_symbols,
    	var_group_positions=var_group_positions,
    	var_group_labels=var_group_labels,
    	standard_scale=standard_scale,
    	var_group_rotation=var_group_rotation,
    	layer=layer,
    	stripplot=stripplot,
    	jitter=jitter,
    	size=size,
    	scale=scale,
    	order=order,
    	swap_axes=swap_axes,
    	show=show,
    	save=save,
    	row_palette=row_palette)
