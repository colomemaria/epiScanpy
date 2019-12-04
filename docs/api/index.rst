.. automodule:: episcanpy

API
===

Import epiScanpy's high-level API as::

   import episcanpy.api as epi

Count Matrices: CT
------------------

Loading data, loading annotations, building count matrices, filtering of lowly covered methylation variables.
Filtering of lowly covered cells.

Building count matrices
~~~~~~~~~~~~~~~~~~~~~~~

Quickly build a count matrix from tsv/tbi file.

.. autosummary::
   :toctree: .

   ct.bld_mtx_fly

Load features
~~~~~~~~~~~~~

In order to build a count matrix for either methylation or open chromatin data, loading the segmentation of the genome of interest or the set of features of interest is a prerequirement.

.. autosummary::
   :toctree: .

   ct.load_features
   ct.make_windows
   ct.size_feature_norm
   ct.plot_size_features
   ct.name_features

Reading methylation file
~~~~~~~~~~~~~~~~~~~~~~~~

Functions to read methylation files, extract methylation and buildthe count matrices:
 
 
.. autosummary::
   :toctree: .

   ct.build_count_mtx
   ct.read_cyt_summary
   ct.load_met_noimput

   
Reading open chromatin(ATAC) file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ATAC-seq specific functions to build count matrices and load data:
 
 
.. autosummary::
   :toctree: .
   
   ct.bld_atac_mtx
   ct.save_sparse_mtx

   
General functions
~~~~~~~~~~~~~~~~~

Functions non -omic specific:

.. autosummary::
   :toctree: .

   ct.save_sparse_mtx



Preprocessing: PP
------------------

Imputing missing data (methylation), filtering lowly covered cells or variables, correction for batch effect.


.. autosummary::
   :toctree: .

   pp.coverage_cells
   pp.commoness_features
   pp.select_var_feature
   pp.binarize
   pp.lazy
   pp.load_metadata
   pp.read_ATAC_10x
   pp.filter_cells
   pp.filter_features
   pp.normalize_total
   pp.pca
   pp.normalize_per_cell
   pp.regress_out
   pp.subsample
   pp.downsample_counts
   pp.neighbors
   pp.sparse



Methylation matrices
~~~~~~~~~~~~~~~~~~~~

Methylation specific count matrices. 

.. autosummary::
   :toctree: .

   pp.imputation_met
   pp.load_met_noimput
   pp.readandimputematrix


   
Tools: TL
---------

.. autosummary::
   :toctree: .

   tl.rank_features
   tl.silhouette
   tl.lazy
   tl.load_markers
   tl.identify_cluster
   tl.top_feature_genes
   tl.find_genes
   tl.diffmap
   tl.draw_graph
   tl.tsne
   tl.umap
   tl.louvain
   tl.leiden



Plotting: PL
------------

The plotting module :class:`episcanpy.plotting` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. autosummary::
   :toctree: .
   

   pl.pca
   pl.pca_floadings
   pl.pca_overview
   pl.pca_variance_ratio

   pl.tsne
   pl.umap
   pl.diffmappl.draw_graph

   pl.rank_feat_groups
   pl.rank_feat_groups_violin
   pl.rank_feat_groups_dotplot
   pl.rank_feat_groups_stacked_violin
   pl.rank_feat_groups_matrixplot
   pl.rank_feat_groups_heatmap
   pl.rank_feat_groups_tracksplot
   pl.cal_var

   pl.violin
   pl.scatter
   pl.ranking
   pl.clustermap
   pl.stacked_violin
   pl.heatmap
   pl.dotplot
   pl.matrixplot
   pl.tracksplot
   pl.dendrogram
   pl.correlation_matrix

   pl.prct_overlap
   pl.overlap_heatmap
   pl.cluster_composition
   pl.silhouette
   pl.silhouette_tot



