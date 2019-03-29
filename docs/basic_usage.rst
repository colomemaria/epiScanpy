Usage Principles
----------------


Import the epiScanpy API as::

    import episcanpy.api as epi
    import anndata as ad

Workflow
^^^^^^^^

First we need to build the count matrix. It requires (sometimes) -omic specific functions.
All the functions to build the count matrices (that are either for ATAC, methylation or other) will  use ``epi.ct``

When it comes to calculationg tSNE, UMAP, PCA etc. we take advantage of the shared datastructure with scanpy and we can use most (if not all) Scanpy functions

To see Scanpy usage principles: <https://scanpy.readthedocs.io/en/latest/basic_usage.html>`__.


The typical workflow consists of subsequent calls of ATAC specific processing tools
in ``epi.ch``, e.g.::

    epi.ch.load_features(file_features, **tool_params)  # to load annotation files (bad example as it works the same for mt

where ``adata`` is an :class:`~anndata.AnnData` object. Each of these calls adds annotation to an expression matrix *X*, which stores *n_obs* observations (cells) of *n_vars* variables (genes). For each tool, there typically is an associated plotting function in ``sc.pl``::

    sc.pl.tsne(adata, **plotting_params)
    
And, just because I can, I will inserate a "random" figure. Truth being told,I would like to put the final pipeline.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"><img src="https://github.com/DaneseAnna/Episcanpy/tree/master/docs/api/umapSatb2_CLUSTER_NORM.png" style="width: 100px"><img src="https://github.com/DaneseAnna/Episcanpy/tree/master/docs/api/umapSatb2_CLUSTER_NORM.png" style="width: 100px"><img src="https://github.com/DaneseAnna/Episcanpy/tree/master/docs/api/umapSatb2_CLUSTER_NORM.png" style="width: 100px"><img src="https://github.com/DaneseAnna/Episcanpy/tree/master/docs/api/umapSatb2_CLUSTER_NORM.png" style="width: 200px">




.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/

