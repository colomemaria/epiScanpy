.. automodule:: episcanpy.api

API
===


Import Scanpy's high-level API as::

   import episcanpy.api as epi

Preprocessing: PP
------------------

Filtering of highly-variable genes, batch-effect correction, per-cell normalization, preprocessing recipes.

Basic Preprocessing
~~~~~~~~~~~~~~~~~~~

For visual quality control, see  :func:`~episcanpy.api.new_basic_functions.load_features` :func:`~episcanpy.api.load.extract_CG` :func:`~episcanpy.api.load.extract_CH` and
 in the :doc:`plotting API <plotting>`.

.. autosummary::
   :toctree: .

   functions.extract_meth.extract_CG