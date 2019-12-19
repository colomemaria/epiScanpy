import numpy as np
import anndata as ad
import pandas as pd
import warnings
from warnings import warn
from scipy.sparse  import  csc_matrix

def load_peak_matrix(matrix_file, path=''):
    
    """
    Deprecated - Use load_atac_matrix instead

    Load existing peak x cell matrix into an anndata object.
    Input existing peak matrix in the format as in Cusanovich et al. 2018 
    First row of the matrix file is considered the header (chr start end annot cell1...cellN)
    Paramters
    ---------
    matrix_file: txt or tsv file with the matrix in shape peaks x cells  
    
    Return
    ------
    AnnData  
    """
    warn('Deprecated - Use load_atac_matrix instead')
    peak_name = []
    cell_matrix = [] 
    with open(path+matrix_file) as f:
        head = f.readline().rstrip("\n").split('\t')
        for line in f:
            line = line.rstrip("\n").split('\t')
            peak_name.append(line[3]) # for some reason it has rownames
            line2 = []
            for x in line[5:]:
                line2.append(int(x))
            cell_matrix.append(line2)
    cell_names = head[4:]
    
    cell_matrix=csc_matrix(np.matrix(cell_matrix)).transpose()   
    adata = ad.AnnData(cell_matrix,
                   obs=pd.DataFrame(index=cell_names),
                   var=pd.DataFrame(index=peak_name))
    adata.uns['omic'] = 'ATAC'
    
    return(adata)

def load_bedtool_matrix(matrix_file, path=''):
    
    """
    Deprecated - Use load_atac_matrix instead

    Load existing peak x cell matrix into an anndata object.
    Input existing peak matrix in the format as in Cusanovich et al. 2018 
    First row of the matrix file is considered the header (chr start end annot cell1...cellN)
    Paramters
    ---------
    matrix_file: txt or tsv file with the matrix in shape peaks x cells  
    
    Return
    ------
    AnnData  
    """
    warn('Deprecated - Use load_atac_matrix instead')
    peak_name = []
    cell_matrix = [] 
    with open(path+matrix_file) as f:
        head = f.readline().rstrip("\n").split('\t')
        for line in f:
            line = line.rstrip("\n").split('\t')
            peak_name.append(line[3]) # for some reason it has rownames
            line2 = []
            for x in line[5:]:
                line2.append(int(x))
            cell_matrix.append(line2)
    cell_names = head[4:]
    
    cell_matrix=csc_matrix(np.matrix(cell_matrix)).transpose()   
    adata = ad.AnnData(cell_matrix,
                   obs=pd.DataFrame(index=cell_names),
                   var=pd.DataFrame(index=peak_name))
    adata.uns['omic'] = 'ATAC'
    
    return(adata)

def load_atac_matrix(matrix_file, path='', compression=None):
    """

    Load existing peak x cell matrix into an anndata object.
    Input existing peak matrix in the format as in Cusanovich et al. 2018 
    First row of the matrix file is considered the header (chr start end annot cell1...cellN)
    Paramters
    ---------
    matrix_file: txt or tsv file with the matrix in shape peaks x cells  
    
    compression: if compression is 'gzip'
    Return
    ------
    AnnData
    """
    
    if compression=='gzip':
        data = pd.read_csv(path+matrix_file,
            sep='\t',
            header=0,
            compression='gzip')
    else:
        data = pd.read_csv(path+matrix_file,
            sep='\t',
            header=0)

    features = data.annot.tolist()
    barcodes = data.columns.tolist()[4:]
    data = csc_matrix(data.iloc[:,4:].values)
    adata = ad.AnnData(data.transpose(),
                   obs=pd.DataFrame(index=barcodes),
                   var=pd.DataFrame(index=features))
    adata.uns['omic'] = 'ATAC'
    return(adata)


