.. automodule:: episcanpy

API
===

Import Scanpy's high-level API as::

   import episcanpy.api as epi

Count Matrices: CT
------------------

Loading data, loading annotations, building count matrices, filtering of lowly covered methylation variables.
Filtering of lowly covered cells.

Building Count matrices
~~~~~~~~~~~~~~~~~~~~~~~

For visual quality control, see  :func:`~episcanpy.pp.extract_CG` :func:`~episcanpy.ct.read_methylation_file` :func:`~episcanpy.ct.make_windows` and  in :mod:`episcanpy.plotting`.
 in the :doc:`plotting API <plotting>`.
 
 
.. autosummary::
   :toctree: .

   ct.read_methylation_file
   ct.load_features
   ct.make_windows
   

Preprocessing: PP
------------------

Imputing missing data (methylation), filterinf lowly covered cells or variables, correction for batch effect...

Preprocessing
~~~~~~~~~~~~~

For the sake of testing the readthedocs, see  :func:`~episcanpy.api.load.load_features` :func:`~episcanpy.load_features.load_features` :func:`~episcanpy.api.load.extract_CH` and  in :mod:`episcanpy.plotting`.
 in the :doc:`plotting API <plotting>`.
 
 
.. autosummary::
   :toctree: .

   pp.extract_CH
   pp.extract_CG
   


Plotting: PL
------------

The plotting module :class:`episcanpy.plotting` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. toctree::
   :hidden:
   :maxdepth: 1

   whatever



