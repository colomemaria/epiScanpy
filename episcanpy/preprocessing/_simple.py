import scanpy as sc
import anndata as ad
from scipy.sparse import csc_matrix

def sparse(adata, copy=False):
	"""
	Transform adata.X from a matrix or array to a csc sparse matrix.
	"""
	if copy:
		adata2 = adata.copy()
		adata2.X = csc_matrix(adata2.X)
		return(adata2)
	else:
		adata.X = csc_matrix(adata.X)

