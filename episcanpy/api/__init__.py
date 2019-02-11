from anndata import AnnData

from anndata import read as read_h5ad
from anndata import read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools

from .. import __version__

#from . import load
#from . import new_basic_function
#from . import read_meth_file
# from ..readwrite import read, read_10x_h5, read_10x_mtx, write, read_params, write_params
# from . import datasets
# from . import export_to
# from . import logging
# from . import queries

from . import pp
from . import ct


# unfortunately, we cannot put this here as long as we have simple global
# variables in settings... they couldn't be set in this case...
# the main drawback is that we have to import set_figure_params
# to show in the docs for that reason...
# it would be nice to make the simple data types "properties of the
# module"... putting setters and getters for all of them wouldn't be very nice
from .. import settings
# for now - or maybe as the permanently favored solution - put the single function here
from ..settings import set_figure_params

# some stuff that is not actually documented...
from .. import utils

import sys
utils.annotate_doc_types(sys.modules[__name__], 'scanpy')
del sys


__doc__ = """\
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

For visual quality control, see  :func:`~episcanpy.pp.extract_CG` :func:`~episcanpy.load.read_methylation_file` :func:`~episcanpy.ct.make_windows` and  in :mod:`episcanpy.plotting`.
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
"""
