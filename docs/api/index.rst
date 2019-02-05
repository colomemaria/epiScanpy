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

For visual quality control, see  :func:`~episcanpy.api.load.load_features` :func:`~episcanpy.load_features.load_features` :func:`~episcanpy.api.load.extract_CH` and  in :mod:`episcanpy.plotting`.
 in the :doc:`plotting API <plotting>`.
 
 
.. autosummary::
   :toctree: .

   load.extract_CH
   load.extract_CG
   load.load_features



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
