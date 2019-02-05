.. automodule:: episcanpy

API
===


Import Scanpy's high-level API as::

   import episcanpy.api as epi

Preprocessing: PP
------------------

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For visual quality control, see  :func:`~episcanpy.api.new_basic_functions.load_features` :func:`~episcanpy.load_features.load_features` :func:`~episcanpy.api.load.extract_CH` and
 in the :doc:`plotting API <plotting>`.

.. autosummary::
   :toctree: .

   functions.extract_meth.extract_CG
   load.extract_CH



Plotting: PL
------------

The plotting module :class:`episcanpy.plotting` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. toctree::
   :hidden:
   :maxdepth: 1

   plotting




Functions: FUN
--------------

The plotting module :class:`episcanpy.functions` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. toctree::
   :hidden:
   :maxdepth: 1

   functions




Reading
-------

*Note:* For reading annotation use :ref:`pandas.read_… <pandas:/io.rst#io-tools-text-csv-hdf5>`
and add it to your :class:`anndata.AnnData` object.
The following read functions are intended for the numeric data in the data matrix `X`.

Read common file formats using

.. autosummary::
   :toctree: .

   read

Read 10x formatted hdf5 files and directories containing `.mtx` files using

.. autosummary::
   :toctree: .

    read_10x_h5
    read_10x_mtx

Read other formats using functions borrowed from :mod:`anndata`

.. autosummary::
   :toctree: .

   read_h5ad
   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools
