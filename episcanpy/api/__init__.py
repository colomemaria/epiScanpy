from anndata import AnnData

from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

from .. import __version__

from . import pp
from . import tl
#from . import functions
from . import ct
from . import pl

MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']

__doc__ = """\
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
   ct.size_feature_norm
   ct.plot_size_features
   ct.name_features

Reading methylation file
~~~~~~~~~~~~~~~~~~~~~~~~

Functions to read methylation files, extract methylation and buildthe count matrices:
 
 
.. autosummary::
   :toctree: .

   ct.build_count_mtx
   ct.read_cyt_summary
   ct.load_met_noimput

   
Reading open chromatin(ATAC) file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ATAC-seq specific functions to build count matrices and load data:
 
 
.. autosummary::
   :toctree: .
   
   ct.bld_atac_mtx
   ct.save_sparse_mtx

   
General functions
~~~~~~~~~~~~~~~~~

Functions non -omic specific:

   .. autosummary::
   :toctree: .

   ct.save_sparse_mtx


Preprocessing: PP
------------------

Imputing missing data (methylation), filtering lowly covered cells or variables, correction for batch effect.


.. autosummary::
   :toctree: .

   pp.coverage_cells
   pp.commoness_features
   pp.binarize
   pp.lazy
   pp.load_metadata
   pp.read_ATAC_10x


Methylation matrices
~~~~~~~~~~~~~~~~~~~~

Methylation specific count matrices. 

.. autosummary::
   :toctree: .

   pp.imputation_met
   pp.load_met_noimput
   pp.readandimputematrix


   
Tools: TL
---------

.. autosummary::
   :toctree: .

   tl.rank_features
   tl.silhouette
   tl.lazy
   tl.load_markers
   tl.identify_cluster


Plotting: PL
------------

The plotting module :class:`episcanpy.plotting` largely parallels the ``tl.*`` and a few of the ``pp.*`` functions.
For most tools and for some preprocessing functions, you'll find a plotting function with the same name.

.. autosummary::
   :toctree: .
   pl.prct_overlap
   pl.overlap_heatmap
   pl.silhouette
   pl.silhouette_tot


   
"""
