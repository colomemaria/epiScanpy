import numpy as np
import anndata as ad
import pandas as pd

def load_peak_matrix(matrix_file, path=''):
    
    """
    Load existing peak x cell matrix into an anndata object.
    Input existing peak matrix in the format as in Cusanovich et al. 
    First row of the matrix file is considered the header (chr start end annot cell1...cellN)
    Paramters
    ---------
    matrix_file: txt or tsv file with the matrix in shape peaks x cells  
    
    Return
    ------
    AnnData
    """
    
    peak_name = []
    cell_matrix = [] 
    with open(path+matrix_file) as f:
        head = f.readline().split('\t')
        head[len(head)-1] = head[len(head)-1].split("\n")[0]
        for line in f:
            line = line.split('\t')
            line[len(line)-1] = line[len(line)-1].split("\n")[0]
            peak_name.append(line[3]) # for some reason it has rownames
            line2 = []
            for x in line[4:]:
                if x =='0':
                    line2.append(0)
                else:
                    line2.append(1)
            cell_matrix.append(line2)
    cell_names = head[4:]
    
    cell_matrix=np.matrix(cell_matrix)
    cell_matrix = cell_matrix.transpose()
    
    adata = ad.AnnData(cell_matrix,
                   obs=pd.DataFrame(index=cell_names),
                   var=pd.DataFrame(index=peak_name))
    adata.uns['omic'] = 'ATAC'
    
    return(adata)