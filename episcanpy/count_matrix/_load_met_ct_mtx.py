import numpy as np
import anndata as ad

def load_met_noimput(matrix_file, path='', save=False):
    """
    read the raw count matrix and convert it into an AnnData object.
    write down the matrix as .h5ad (AnnData object) if save = True.
    Return AnnData object
    """
    matrix = []
    cell_names = []
    feature_names = []
    with open('../../enhancer_c1_CG_paper_unique_corrected2.txt') as f:
        line = f.readline()[:-2].split('\t')
        if line[0] == 'sample_name':
            feature_names = line[1:-1]
        else:
            matrix.append(line[1:])
            cell_names.append(line[0])
        if matrix == []:
            line = f.readline()[:-2].split('\t')
            matrix.append(line[1:])
            cell_names.append(line[0])
        for line in f:
            line = line[:-2].split('\t')
            matrix.append(line[1:])
            cell_names.append(line[0])

    matrix = np.array(matrix)
    
    if feature_names != []:
        adata = ad.AnnData(matrix, obs=pd.DataFrame(index=cell_names),  var=pd.DataFrame(index=feature_names))
    else:
        adata = ad.AnnData(matrix, obs=pd.DataFrame(index=cell_names))
    
    adata.uns['omic'] = 'methylation'
    
    if save:
        adata.write("".join([".".split(matrix_file)[0],'.h5ad']))
        
    return(adata)
