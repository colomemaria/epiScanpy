import numpy as np
from scipy.io import mmread
import pandas as pd
import anndata as ad


def read_ATAC_10x(matrix, cell_names='', var_names='', path_file=''):
    """
    Load sparse matrix (including matrices corresponding to 10x data) as AnnData objects.
    read the mtx file, tsv file coresponding to cell_names and the bed file containing the variable names

    Parameters
    ----------
    matrix: sparse count matrix

    cell_names: optional, tsv file containing cell names

    var_names: optional, bed file containing the feature names

    Return
    ------
    AnnData object

    """

    
    mat = mmread(''.join([path_file, matrix]))
    mat = mat.toarray()
    mat = np.matrix(mat.transpose())
    
    with open(path_file+cell_names) as f:
        barcodes = f.readlines() 
        barcodes = [x[:-1] for x in barcodes]
        
    with open(path_file+var_names) as f:
        var_names = f.readlines()
        var_names = ["_".join(x[:-1].split('\t')) for x in var_names]
        
    adata = ad.AnnData(mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=var_names))
    adata.uns['omic'] = 'ATAC'
    
    return(adata)
