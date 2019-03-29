import numpy as np
from scipy.io import mmread
import pandas as pd
import anndata as ad


def read_ATAC(matrix, cell_names='', var_names='', path_file=''):
    # read the mtx file, the tsv file cell_names and the bed file for var_names"
    
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
