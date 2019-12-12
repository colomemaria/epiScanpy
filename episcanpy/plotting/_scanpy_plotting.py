## everything here is copied/inspired from scanpy

# to add: violin, ranking, clustermap, stacked_violin, heatmap, dotplot,
# matrixplot, tracksplot, dendrogram, correlation_matrix
import scanpy as sc
from anndata import AnnData
from typing import Optional, Union  # Special
from typing import Sequence, Collection, Iterable  # ABCs
from typing import Tuple, List  # Classes
from matplotlib.axes import Axes
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib import patheffects
from matplotlib.colors import is_color_like, Colormap, ListedColormap


def pca(adata,
	color=None,
	feature_symbols=None,
	use_raw=None,
	layer=None,
	sort_order=True,
	groups=None,
	components=None,
	projection='2d',
	legend_loc='right margin',
	legend_fontsize=None,
	legend_fontweight=None,
	legend_fontoutline=None,
	size=None,
	color_map=None,
	palette=None,
	frameon=None,
	ncols=None,
	wspace=None,
	hspace=0.25,
	title=None,
	return_fig=None,
	show=None,
	save=None):
	"""\
    Scatter plot in PCA coordinates.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """

	sc.pl.pca(adata, color=color, gene_symbols=feature_symbols, use_raw=use_raw,
        layer=layer, sort_order=sort_order, groups=groups, components=components,
        projection=projection, legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight, legend_fontoutline=legend_fontoutline,
        size=size, color_map=color_map, palette=palette, frameon=frameon, ncols=ncols,
        wspace=wspace, hspace=hspace, title=title, return_fig=return_fig, show=show,
        save=save)

def tsne(adata,
	color=None,
	feature_symbols=None,
	use_raw=None,
	layer=None,
	sort_order=True,
	groups=None,
	components=None,
	projection='2d',
	legend_loc='right margin',
	legend_fontsize=None,
	legend_fontweight=None,
	legend_fontoutline=None,
	size=None,
	color_map=None,
	palette=None,
	frameon=None,
	ncols=None,
	wspace=None,
	hspace=0.25,
	title=None,
	return_fig=None,
	show=None,
	save=None):
	"""\
    Scatter plot in tSNE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
	sc.pl.tsne(adata, color=color, gene_symbols=feature_symbols, use_raw=use_raw,
        layer=layer, sort_order=sort_order, groups=groups, components=components,
        projection=projection, legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight, legend_fontoutline=legend_fontoutline,
        size=size, color_map=color_map, palette=palette, frameon=frameon, ncols=ncols,
        wspace=wspace, hspace=hspace, title=title, return_fig=return_fig, show=show,
        save=save)

def umap(adata,
	color=None,
	feature_symbols=None,
	use_raw=None,
	layer=None,
	sort_order=True,
	groups=None,
	components=None,
	projection='2d',
	legend_loc='right margin',
	legend_fontsize=None,
	legend_fontweight=None,
	legend_fontoutline=None,
	size=None,
	color_map=None,
	palette=None,
	frameon=None,
	ncols=None,
	wspace=None,
	hspace=0.25,
	title=None,
	return_fig=None,
	show=None,
	save=None):
	"""\
    Scatter plot in UMAP basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
	sc.pl.umap(adata, color=color, gene_symbols=feature_symbols, use_raw=use_raw,
        layer=layer, sort_order=sort_order, groups=groups, components=components,
        projection=projection, legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight, legend_fontoutline=legend_fontoutline,
        size=size, color_map=color_map, palette=palette, frameon=frameon, ncols=ncols,
        wspace=wspace, hspace=hspace, title=title, return_fig=return_fig, show=show,
        save=save)


def diffmap(adata,
	color=None,
	feature_symbols=None,
	use_raw=None,
	layer=None,
	sort_order=True,
	groups=None,
	components=None,
	projection='2d',
	legend_loc='right margin',
	legend_fontsize=None,
	legend_fontweight=None,
	legend_fontoutline=None,
	size=None,
	color_map=None,
	palette=None,
	frameon=None,
	ncols=None,
	wspace=None,
	hspace=0.25,
	title=None,
	return_fig=None,
	show=None,
	save=None):
	"""\
    Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
	sc.pl.diffmap(adata, color=color, gene_symbols=feature_symbols, use_raw=use_raw,
        layer=layer, sort_order=sort_order, groups=groups, components=components,
        projection=projection, legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight, legend_fontoutline=legend_fontoutline,
        size=size, color_map=color_map, palette=palette, frameon=frameon, ncols=ncols,
        wspace=wspace, hspace=hspace, title=title, return_fig=return_fig, show=show,
        save=save)

def draw_graph(adata,
    layout=None,
	color=None,
	feature_symbols=None,
	use_raw=None,
	layer=None,
	sort_order=True,
	groups=None,
	components=None,
	projection='2d',
	legend_loc='right margin',
	legend_fontsize=None,
	legend_fontweight=None,
	legend_fontoutline=None,
	size=None,
	color_map=None,
	palette=None,
	frameon=None,
	ncols=None,
	wspace=None,
	hspace=0.25,
	title=None,
	return_fig=None,
	show=None,
	save=None):
	"""\
    Scatter plot in graph-drawing basis.

    Parameters
    ----------
    {adata_color_etc}
    layout : {{`'fa'`, `'fr'`, `'drl'`, ...}}, optional (default: last computed)
        One of the :func:`~scanpy.tl.draw_graph` layouts.
        By default, the last computed layout is used.
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.
    """
	sc.pl.draw_graph(adata, layout, color=color, gene_symbols=feature_symbols, use_raw=use_raw,
        layer=layer, sort_order=sort_order, groups=groups, components=components,
        projection=projection, legend_loc=legend_loc, legend_fontsize=legend_fontsize,
        legend_fontweight=legend_fontweight, legend_fontoutline=legend_fontoutline,
        size=size, color_map=color_map, palette=palette, frameon=frameon, ncols=ncols,
        wspace=wspace, hspace=hspace, title=title, return_fig=return_fig, show=show,
        save=save)

def rank_feat_groups_violin(adata,
    groups=None,
    n_features=20,
    feature_names=None,
    feature_symbols=None,
    use_raw=None,
    key='rank_features_groups',
    split=True,
    scale='width',
    strip=True,
    jitter=True,
    size=1,
    ax=None,
    show=None,
    save=None):
	"""\
    Plot ranking of features for all tested comparisons.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        List of group names.
    n_features
        Number of features to show. Is ignored if `feature_names` is passed.
    feature_names
        List of features to plot. Is only useful if interested in a custom feature list,
        which is not the result of :func:`epi.tl.rank_features`.
    feature_symbols
        Key for field in `.var` that stores feature symbols if you do not want to
        use `.var_names` displayed in the plot.
    use_raw : `bool`, optional (default: `None`)
        Use `raw` attribute of `adata` if present. Defaults to the value that
        was used in :func:`~scanpy.tl.rank_genes_groups`.
    split
        Whether to split the violins or not.
    scale
        See :func:`~seaborn.violinplot`.
    strip
        Show a strip plot on top of the violin plot.
    jitter
        If set to 0, no points are drawn. See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    {show_save_ax}
    """
	sc.pl.rank_genes_groups_violin(adata, groups, n_genes=n_features,
		gene_names=feature_names, gene_symbols=feature_symbols, use_raw=use_raw,
		key=key, split=split, scale=scale, strip=strip, jitter=jitter, size=size, 
        ax=ax, show=show, save=save)

def rank_feat_groups(adata,
    groups = None,
    n_features= 20,
    feature_symbols = None,
    key = 'rank_features_groups',
    fontsize = 8,
    ncols= 4,
    sharey= True,
    show = None,
    save = None,
    ax = None,):
	"""\
    Plot ranking of features.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the feature ranking.
    feature_symbols
        Key for field in `.var` that stores feature symbols if you do not want to
        use `.var_names`.
    n_features
        Number of feature to show.
    fontsize
        Fontsize for feature names.
    ncols
        Number of panels shown per row.
    sharey
        Controls if the y-axis of each panels should be shared. But passing
        `sharey=False`, each panel has its own y-axis range.
    {show_save_ax}
    """
	sc.pl.rank_genes_groups(adata, groups, n_genes=n_features, gene_symbols=feature_symbols,
		key=key, fontsize=fontsize, ncols=ncols, sharey=sharey, show=show, save=save , ax=ax)

def rank_feat_groups_dotplot(adata,
    groups = None,
    n_features = 10,
    groupby = None,
    key ='rank_features_groups',
    show = None,
    save = None):
	"""\
    Plot ranking of features using dotplot plot (see :func:`~scanpy.pl.dotplot`)

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the feature ranking.
    n_features
        Number of features to show.
    groupby
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank features groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.dotplot`).
    key
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds
        Are passed to :func:`~scanpy.pl.dotplot`.
    """
	sc.pl.rank_genes_groups_dotplot(adata, groups, n_genes=n_features, groupby=groupby,
        key=key, show=show, save=save)

def rank_feat_groups_heatmap(adata,
    groups= None,
    n_features = 10,
    groupby = None,
    key = 'rank_features_groups',
    show = None,
    save = None,):
	"""\
    Plot ranking of features using heatmap plot (see :func:`~scanpy.pl.heatmap`)

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the feature ranking.
    n_features
        Number of features to show.
    groupby
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank features groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.heatmap`).
    key
        Key used to store the ranking results in `adata.uns`.
    **kwds
        Are passed to :func:`~scanpy.pl.heatmap`.
    {show_save_ax}
    """
	sc.pl.rank_genes_groups_heatmap(adata, groups, n_genes=n_features, groupby=groupby,
        key=key, show=show, save=save)

def rank_feat_groups_stacked_violin(adata,
    groups = None,
    n_features = 10,
    groupby = None,
    key = 'rank_features_groups',
    show = None,
    save = None):
	"""\
    Plot ranking of features using stacked_violin plot (see :func:`~scanpy.pl.stacked_violin`)

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups : `str` or `list` of `str`
        The groups for which to show the feature ranking.
    n_features : `int`, optional (default: 10)
        Number of features to show.
    groupby : `str` or `None`, optional (default: `None`)
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank features groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.stacked_violin`).
    key
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds
        Are passed to :func:`~scanpy.pl.stacked_violin`.
    """
	sc.pl.rank_genes_groups_stacked_violin(adata, groups, n_genes=n_features, groupby=groupby,
        key=key, show=show, save=save)

def rank_feat_groups_matrixplot(adata,
    groups = None,
    n_features = 10,
    groupby = None,
    key = 'rank_features_groups',
    show = None,
    save = None):
	"""\
    Plot ranking of features using matrixplot plot (see :func:`~scanpy.pl.matrixplot`)

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the feature ranking.
    n_features
        Number of features to show.
    groupby
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank features groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.matrixplot`).
    key
        Key used to store the ranking results in `adata.uns`.
    {show_save_ax}
    **kwds
        Are passed to :func:`~scanpy.pl.matrixplot`.
    """
	sc.pl.rank_genes_groups_matrixplot(adata, groups, n_genes=n_features, groupby=groupby,
		key=key, show=show, save=save)

def rank_feat_groups_tracksplot(adata,
    groups = None,
    n_features = 10,
    groupby = None,
    key = 'rank_features_groups',
    show = None,
    save = None):
	"""\
    Plot ranking of features using heatmap plot (see :func:`~scanpy.pl.heatmap`)

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        The groups for which to show the feature ranking.
    n_features
        Number of features to show.
    groupby
        The key of the observation grouping to consider. By default,
        the groupby is chosen from the rank features groups parameter but
        other groupby options can be used.  It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories` (see :func:`~scanpy.pl.heatmap`).
    key
        Key used to store the ranking results in `adata.uns`.
    **kwds
        Are passed to :func:`~scanpy.pl.tracksplot`.
    {show_save_ax}
    """
	sc.pl.rank_genes_groups_tracksplot(adata, groups, n_genes=n_features, groupby=groupby,
        key=key, show=show, save=save)

def pca_loadings(adata,
    components = None,
    include_lowest = True,
    show = None,
    save = None):
	"""\
    Rank features according to contributions to PCs.

    Parameters
    ----------
    adata
        Annotated data matrix.
    components
        For example, ``'1,2,3'`` means ``[1, 2, 3]``, first, second, third
        principal component.
    include_lowest
        Show the features with both highest and lowest loadings.
    show
        Show the plot, do not return axis.
    save
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    """
	sc.pl.pca_loadings(adata, components, include_lowest, show, save)


def pca_overview(adata, color=None, use_raw=True, show=None, save=None):
	"""\
    Plot PCA results.

    The parameters are the ones of the scatter plot. Call pca_ranking separately
    if you want to change the default settings.

    Parameters
    ----------
    
    adata
        Annotated data matrix.
    color : string or list of strings, optional (default: `None`)
        Keys for observation/cell annotation either as list `["ann1", "ann2"]` or
        string `"ann1,ann2,..."`.
    use_raw : `bool`, optional (default: `True`)
        Use `raw` attribute of `adata` if present.
    {scatter_bulk}
    show : bool, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {{`'.pdf'`, `'.png'`, `'.svg'`}}.
    """
	sc.pl.pca_overview(adata, color=color, use_raw=use_raw, show=show, save=save)

def pca_variance_ratio(adata, n_pcs=30, log=False, show=None, save=None):
	"""\
    Plot the variance ratio.

    Parameters
    ----------
    n_pcs : `int`, optional (default: `30`)
         Number of PCs to show.
    log : `bool`, optional (default: `False`)
         Plot on logarithmic scale..
    show : `bool`, optional (default: `None`)
         Show the plot, do not return axis.
    save : `bool` or `str`, optional (default: `None`)
        If `True` or a `str`, save the figure.
        A string is appended to the default filename.
        Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    """
	sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=log, show=show, save=save)



def correlation_matrix(
    adata: AnnData,
    groupby: str,
    show_correlation_numbers: bool = False,
    dendrogram: Union[bool, str, None] = None,
    figsize: Optional[Tuple[float, float]] = None,
    show: Optional[bool] = None,
    save: Optional[Union[bool, str]] = None,
    ax: Optional[Axes] = None,
    **kwds,
) -> Union[Axes, List[Axes]]:
    """Plots the correlation matrix computed as part of `sc.tl.dendrogram`.

    Parameters
    ----------
    adata
    groupby
        Categorical data column used to create the dendrogram
    show_correlation_numbers
        If `show_correlation` is True, plot the correlation number on top of each cell.
    dendrogram
        If True or a valid dendrogram key, a dendrogram based on the hierarchical clustering
        between the `groupby` categories is added. The dendrogram information is computed
        using :func:`scanpy.tl.dendrogram`. If `tl.dendrogram` has not been called previously
        the function is called with default parameters.
    figsize
        By default a figure size that aims to produce a squared correlation matrix plot is used.
        Format is (width, height)
    {show_save_ax}
    **kwds
        Only if `show_correlation` is True:
        Are passed to :func:`matplotlib.pyplot.pcolormesh` when plotting the
        correlation heatmap. Useful values to pas are `vmax`, `vmin` and `cmap`.

    Returns
    -------

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, 'bulk_labels')
    >>> sc.pl.correlation(adata, 'bulk_labels')
    """
    sc.pl.correlation_matrix(
    adata=adata,
    groupby=groupby,
    show_correlation_numbers=show_correlation_numbers,
    dendrogram=dendrogram,
    figsize=figsize,
    show=show,
    save=save,
    ax=ax,
    **kwds)

def dendrogram(
    adata: AnnData,
    groupby: str,
    dendrogram_key: Optional[str] = None,
    orientation: str = 'top',
    remove_labels: bool = False,
    show: Optional[bool] = None,
    save: Optional[bool] = None,
):
    """Plots a dendrogram of the categories defined in `groupby`.

    See :func:`~scanpy.tl.dendrogram`.

    Parameters
    ----------
    adata
    groupby
        Categorical data column used to create the dendrogram
    dendrogram_key
        Key under with the dendrogram information was stored.
        By default the dendrogram information is stored under .uns['dendrogram_' + groupby].
    orientation
        Options are `top` (default), `bottom`, `left`, and `right`.
        Only when `show_correlation` is False.
    remove_labels
    {show_save_ax}

    Returns
    -------
    :class:`matplotlib.axes.Axes`

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.dendrogram(adata, 'bulk_labels')
    >>> sc.pl.dendrogram(adata, 'bulk_labels')
    """
    dendrogram(
    adata=adata,
    groupby=groupby,
    dendrogram_key=dendrogram_key,
    orientation=orientation,
    remove_labels=remove_labels,
    show=show,
    save=save)
