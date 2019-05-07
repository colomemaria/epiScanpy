from anndata import AnnData

from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

from .. import __version__

from . import pp
from . import tl
#from . import functions
from . import ct
from . import pl

__doc__ = """\
API
===


Import Scanpy's high-level API as::

   import episcanpy.api as epi

Count Matrices: CT
------------------

Loading data, loading annotations, building count matrices, filtering of lowly covered methylation variables.
Filtering of lowly covered cells.

Load & Generate Features
~~~~~~~~~~~~~~~~~~~~~~~~

To load known annotations of the genome to build count matrices. Or, to generate non-overlapping windows of different sizes.

.. autosummary::
   :toctree: .
   
   ct.load_features
   ct.make_windows
   ct.size_feature_norm
   ct.plot_size_features

Building Count matrices
~~~~~~~~~~~~~~~~~~~~~~~

To build the count matrix (based on annotations or windows)
 
.. autosummary::
   :toctree: .

   ct.build_count_mtx
   ct.read_cyt_summary
   

Preprocessing: PP
------------------

Imputing missing data (methylation), filterinf lowly covered cells or variables, correction for batch effect...

Loading the data
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .
   
   pp.readandimputematrix

Quality Controls
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: .
   
   pp.coverage_cells
   pp.commoness_features


Preprocessing
~~~~~~~~~~~~~
 
 
.. autosummary::
   :toctree: .
   
   pp.binarize
   pp.lazy


Tools: TL
---------

These functions should probably be in preprocessing
.. autosummary::
   :toctree: .

   tl.read_ATAC
   tl.readandimputematrix
   tl.read_MET
   
   
tl.lazy should be there and not in preprocessing

.. autosummary::
   :toctree: .
   tl.lazy
   tl. rank_features
   tl.silhouette

Plotting: PL
------------

Most of the functions that actually produce plot.

I need to put back the coverage and commonness functions from above (or a "plotting onluy' equivalent).

.. toctree::
   :hidden:
   :maxdepth: 1

   pl.plot_rank_features
   pl.overlap_heatmap
   pl.prct_overlap
   pl.silhouette
   pl.silhouette_tot
   
   
"""
