import scanpy as sc
import anndata as ad
from scipy.sparse import csc_matrix, csr_matrix
from warnings import warn

def sparse(adata, sparse_format='csr', copy=False):
	"""
	Transform adata.X from a matrix or array to a sparse matrix.
	sparse_format : accepted input are csr or csc
	if copy is True return a new AnnData object. Else overwrite the input object
	"""
	if copy:
		adata2 = adata.copy()
		if sparse_format=='csc':
			adata2.X = csc_matrix(adata2.X)
		elif sparse_format=='csr':
			adata2.X = csr_matrix(adata2.X)
		else:
			warn(sparse_format, ' is not an accepted format. Please specify csc or csr.')
		return(adata2)
	else:
		if sparse_format=='csc':
			adata.X = csc_matrix(adata.X)
		elif sparse_format=='csr':
			adata.X = csr_matrix(adata.X)
		else:
			warn(sparse_format, ' is not an accepted format. Please specify csc or csr.')


