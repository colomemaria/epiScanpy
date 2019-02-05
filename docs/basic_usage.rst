Usage Principles
----------------

We reuse a lot of the Scanpy functions "out of the box". To this end we need to import Scanpy and Anndata

Import the Scanpy API as::

    import episcanpy.api as epi
    import anndata as ad

Workflow
^^^^^^^^

First we need to build the count matrix. It requires (sometimes) -omic specific functions.
Preprocessing functions are therefore divides into 2 categories: ``epi.ch`` for single cell ATAC data and ``epi.mt`` for the methylation data.

When it comes to calculationg tSNE, UMAP, PCA etc. we take advantage of the shared datastructure with scanpy and we can use most (if not all) Scanpy functions

To see Scanpy usage principles: <https://scanpy.readthedocs.io/en/latest/basic_usage.html>`__.


The typical workflow consists of subsequent calls of ATAC specific processing tools
in ``epi.ch``, e.g.::

    epi.ch.load_features(file_features, **tool_params)  # to load annotation files (bad example as it works the same for mt

where ``adata`` is an :class:`~anndata.AnnData` object. Each of these calls adds annotation to an expression matrix *X*, which stores *n_obs* observations (cells) of *n_vars* variables (genes). For each tool, there typically is an associated plotting function in ``sc.pl``::

    sc.pl.tsne(adata, **plotting_params)


.. _Seaborn: http://seaborn.pydata.org/
.. _matplotlib: http://matplotlib.org/
