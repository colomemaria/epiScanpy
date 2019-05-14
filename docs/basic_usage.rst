Usage Principles
----------------


Import the epiScanpy API as::

    import episcanpy.api as epi
    import anndata as ad

Workflow
^^^^^^^^

First step is to build the count matrix. It requires -omic specific approaches as the data have -omic specific limitations.
All the functions to build the count matrices (that are either for ATAC, methylation or other) will  use ``epi.ct``

When it comes to calculationg tSNE, UMAP, PCA etc. we take advantage of the shared datastructure with scanpy and we can use most (if not all) Scanpy functions.

To see Scanpy usage principles: <https://scanpy.readthedocs.io/en/latest/basic_usage.html>`__.


The typical workflow consists of subsequent calls of ATAC specific processing tools
in ``epi.ct``, e.g.::

    epi.ct.load_features(file_features, **tool_params)  # to load annotation files 
    

Second step is to load the matrix and annotations for quality controls, filtering and normalisation. The count matrix and additional informations are stored as ``adata`` is an :class:`~anndata.AnnData` object. 
All functions for quality control and preprocessing are called using ``epi.pp``

To visualise how common features are and what the coverage distribution of the count matrix features, uses: ::
    
    epi.pp.commoness_features(adata, **plotting_params)
    epi.pp.coverage_cells(adata, **plotting_params)
    

To obtain cell-cell distance calculations or low dimensional representation we use the tools developped by  ``scanpy`` and make use of the the ``adata`` object, to store *n_obs* observations (cells) of *n_vars* variables (expression, methylation, chromatin features). For each tool, there typically is an associated plotting function in ``sc.tl`` and``sc.pl``::

        sc.tl.tsne(adata, **tool_params)
        sc.pl.tsne(adata, **plotting_params)
    
Data structure
^^^^^^^^^^^^^^

Similarly to Scanpy, the methylation and ATAC matrices are stored as anndata object. 
    
.. raw:: html

   <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px" style="width: 400px"><img
   <img src="http://github.com/DaneseAnna/episcanpy-pictures/blob/master/heatmap.png" style="width: 300px" style="width: 400px">
  


