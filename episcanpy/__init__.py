# some technical stuff
import sys
from .utils import check_versions, annotate_doc_types
from ._version import get_versions  # version generated by versioneer

__author__ = ', '.join([
    'Anna Danese',
    'Maria Richter',
    'Kridsadakorn Chaichoompu'
])
__email__ = ', '.join([
    'anna.danese@helmholtz-muenchen.de',
    'chaichoompu@helmholtz-muenchen.de',
    # We don’t need all, the main authors are sufficient.
])
__version__ = get_versions()['version']

check_versions()
annotate_doc_types(sys.modules[__name__], 'episcanpy')
del get_versions, sys, check_versions, annotate_doc_types

# the actual API
from ._settings import settings, Verbosity
from . import tools as tl
from . import preprocessing as pp
from . import count_matrix as ct
from . import plotting as pl
# import read functions 
from .preprocessing._load_matrix import read_ATAC_10x, read_h5
from .preprocessing._load_matrix import read_h5 as read_h5_atac

# read functions from AnnData
from anndata import AnnData
from anndata import read
from anndata import read_h5ad, read_csv, read_excel, read_hdf, read_loom, read_mtx, read_text, read_umi_tools




# import multidata information
from ._multidata import MultiData, read_multidata

# import settings for plots 
import scanpy as sc
from sc import set_figure_params

