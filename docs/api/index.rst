.. automodule:: episcanpy

API
===

Import epiScanpy's high-level API as::

   import episcanpy.api as epi

Count Matrices: CT
------------------

Loading data, loading annotations, building count matrices, filtering of lowly covered methylation variables.
Filtering of lowly covered cells.

Load features
~~~~~~~~~~~~~

In order to build a count matrix for either methylation or open chromatin data, loading the segmentation of the genome of interest or the set of features of interest is a prerequirement.

.. autosummary::
   :toctree: .

   ct.load_features
   ct.make_windows

Reading methylation file
~~~~~~~~~~~~~~~~~~~~~~~~

For reading methylation files (ATAC section coming below), extracting methylation and building the count matrix:
(work in progress)
 
 
.. autosummary::
   :toctree: .

   ct.read_methylation_file
   ct.methylation_level
   ct.write_not_sparse_meth
   ct.extract_methylation
   ct.extract_feature_names
   ct.read_meth_fileCG
   ct.read_meth_fileCH
   ct.read_meth_file
   ct.filter_and_average_features
   ct.filter_and_average_features_chrm
   ct.prep_methlevels
   ct.write_list
   ct.write_methlevel
   
Reading open chromatin(ATAC) file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Loading the fucntions and documentation to come. 
(work in progress)
 
 
.. autosummary::
   :toctree: .
   

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



