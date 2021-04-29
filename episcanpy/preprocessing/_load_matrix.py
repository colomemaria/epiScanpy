import h5py
from scipy.sparse import csr_matrix
import pandas as pd
import anndata as ad
import numpy as np
from scipy.io import mmread


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

    
    mat = mmread(''.join([path_file, matrix])).tocsr().transpose()
    
    with open(path_file+cell_names) as f:
        barcodes = f.readlines() 
        barcodes = [x[:-1] for x in barcodes]
        
    with open(path_file+var_names) as f:
        var_names = f.readlines()
        var_names = ["_".join(x[:-1].split('\t')) for x in var_names]
        
    adata = ad.AnnData(mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=var_names))
    adata.uns['omic'] = 'ATAC'
    
    return(adata)


def read_h5(filename, omic='ATAC'):

    """\
    Read `.h5` (hdf5) file.
    
    Parameters
    ----------
    filename
        Filename of data file.
    omic
        Type of single cell data present in the data file

    """
    # open the file
    f = h5py.File(filename, 'r')
    
    #extract obs and var
    barcodes = [x.decode('UTF-8') for x in f['matrix']['barcodes']]
    var = pd.DataFrame()
    if "id" in list(f['matrix']['features'].keys()):
        var["id"] = [x.decode('UTF-8') for x in f['matrix']['features']["id"]]
        var.index = var["id"]
    if "name" in list(f['matrix']['features'].keys()):
        var["name"] = [x.decode('UTF-8') for x in f['matrix']['features']["name"]]
    del var["id"]
    
    
    # extract and re-construct the csr_matrix
    keys = [k for k in f['matrix'].keys()]
    X = f['matrix']["data"][()]
    # try to find row and column names
    rows_cols = [{}, {}]
    for iname, name in enumerate(["indptr", "indices"]):
        if name in keys:
            rows_cols[iname][name] = f['matrix'][name][()]
    
    X = csr_matrix((X, rows_cols[1]["indices"], rows_cols[0]["indptr"]))
    adata = ad.AnnData(X, obs=pd.DataFrame(index=barcodes), var=var)
    adata.uns['omic'] = omic
    return(adata)
    