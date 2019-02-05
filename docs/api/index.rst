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
