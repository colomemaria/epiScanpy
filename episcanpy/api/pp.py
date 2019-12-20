from ..preprocessing._extract import extract_CG, extract_CH
from ..preprocessing._readimpute import load_met_noimput, imputation_met
from ..preprocessing._readimpute import readandimputematrix
from ..preprocessing._quality_control import coverage_cells, commonness_features, binarize
from ..preprocessing._quality_control import correlation_pc, coverage_features, density_features
from ..preprocessing._quality_control import select_var_feature, cal_var, variability_features
from ..preprocessing._recipe import lazy
from ..preprocessing._metadata import load_metadata
from ..preprocessing._load_matrix import read_ATAC_10x
from ..preprocessing._load_atac import load_peak_matrix, load_bedtool_matrix, load_atac_matrix
from ..preprocessing._scanpy_fct import filter_cells, filter_features, normalize_total #, scale
from ..preprocessing._scanpy_fct import pca, normalize_per_cell, regress_out, subsample, downsample_counts 
from ..preprocessing._neighbors import neighbors
from ..preprocessing._simple import sparse
from scanpy.api.pp import log1p