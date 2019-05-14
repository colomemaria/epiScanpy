Usage Principles
----------------


Import the epiScanpy API as::

    import episcanpy.api as epi
    import anndata as ad

Workflow
^^^^^^^^

The first step is to build the count matrix. Because  single-cell epigenomic data types have different characteristics (count data in ATAC-seq versus methylation level in DNA methylation, for example), epiScanpy implements -omic specific approaches to build the count matrix.
All the functions to build the count matrices (for ATAC, methylation or other) will  use ``epi.ct`` (ct = count).

The first step is to load an annotation and then build the count matrix that will be either methylation or ATAC-seq specific. For example using ``epi.ct``, e.g.::

    epi.ct.load_features(file_features, **tool_params)  # to load annotation files 
    epi.ct.build_count_mtx(cell_file_names, omic="ATAC") # to build the ATAC-seq count matrix
    

If you have an already build matrix, you can load it with any additional metadata (such as cell annotations or batches). 

The count matrix, either the one that has been constructed or uploaded, with any additional informations (such as cell annotations or batches) are stored as an :class:`~anndata.AnnData` object. All functions for quality control and preprocessing are called using ``epi.pp`` (pp = preprocessing).

To visualise how common features are and what is the coverage distribution of the count matrix features, use: ::
    
    epi.pp.commoness_features(adata, **processing_params)
    epi.pp.coverage_cells(adata, **processing_params)
    
    
The next step, is the calculation of tSNE, UMAP, PCA etc. For that, we take advantage of the embedding into Scanpy and we use mostly Scanpy functions, which are called using ``sc.tl`` (tl = tool). For that, see Scanpy usage principles: <https://scanpy.readthedocs.io/en/latest/basic_usage.html>`__. For example, to obtain cell-cell distance calculations or low dimensional representation we make use of the ``adata`` object, and store *n_obs* observations (cells) of *n_vars* variables (expression, methylation, chromatin features). For each tool, there typically is an associated plotting function in ``sc.tl`` and ``sc.pl`` (pl = plot) ::

        sc.tl.tsne(adata, **tool_params)
        sc.pl.tsne(adata, **plotting_params)
        
There are also epiScanpy specific tools and plotting functions that can be accessed using ``epi.tl`` and ``epi.pl`` ::

        epi.tl.silhoouette(adata, **tool_params)
        epi.pl.silhouette(adata, **plotting_params)
        epi.pl.prct_overlap(adata, **plotting_params)
        
    
Data structure
^^^^^^^^^^^^^^

Similarly to Scanpy, the methylation and ATAC-seq matrices are stored as anndata object. For more information on the datastructure see here`here <https://anndata.readthedocs.io/en/latest/>`__
    
.. raw:: html

   <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 400px">
  


