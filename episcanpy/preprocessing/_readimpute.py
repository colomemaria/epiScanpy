import numpy as np
import anndata as ad
import pandas as pd

def load_met_noimput(matrix_file, path='', filter_empty=True, save=False):
    """
    read the raw count matrix and convert it into an AnnData object.

    
    Parameters
    ----------
    matrix_file : name of the input file. Required separator is '\t'. 

    path : path to the input file. 
    
    save : if True, write down the matrix as .h5ad (AnnData object). 
    if str specified, save the file with the specified name
   
    Return
    ------
    Return AnnData object

    """
    adata = ad.AnnData(pd.read_csv(path+matrix_file, sep='\t', index_col=0))
    if filter_empty==True:    
        bool_check = []
        for col in range(0,len(adata.var_names)):
            if False not in set(np.isnan(adata.X[:,col].tolist())):
                bool_check.append('remove')
            else:
                bool_check.append('keep')
        adata.var['bool_check'] = bool_check
        adata = adata[:,adata.var['bool_check']!='remove'].copy()
        del adata.var['bool_check']
    
    adata.uns['omic'] = 'methylation'
    adata.uns['imputation'] = 'no_imputation'
    
    if save==True:
        adata.write("".join([".".split(matrix_file)[0],'.h5ad']))
    elif save!=False:
        adata.write(".".join([save.rstrip('.h5ad'), 'h5ad']))

    return(adata)

def imputation_met(adata, number_cell_covered=10, imputation_value='mean', save=None, copy=False):
    """
    Impute missing values in methyaltion level matrices. The imputsation is based on the average
    methylation value of the given variable.
    It also filter out variables that are covered in an unsufficient number of cells in order to 
    reduce the feature space to meaningful variables and discard potential coverage biases. 

    Parameters
    ----------
    adata: AnnData object containing 'nan'

    number_cell_covered: minimum number of cells to be covered in order to retain a variable

    imputation_value: imputation of the missing value can be made either on the mean or the median

    Return
    ------
    Return a new AnnData object

    
    
    """

    # This step need to be sped up and could be multithread.
    # Only the mean available for now. And only the minimum number of cells covered and not the variety of the 
    # methylation levels
    # also, it odes not return the variable annoations and force to add 2 values
    old_features = adata.var_names.tolist()
    
    new_matrix = []
    new_features_name = []
    means = []
    medians = []
    feat_nb = 0

    length1 = len(adata.X[0,:])
    length2 = len(adata.X[:,0])
    adata.obs['coverage_cells'] = [length1 - np.isnan(line).sum() for line in adata.X]
    adata.obs['mean_cell_methylation'] = [np.nansum(line)/length1 for line in adata.X]
    adata.var['coverage_feature'] = [length2 - np.isnan(line).sum() for line in adata.X.T]
    adata.var['mean_feature_methylation'] = [np.nansum(line)/length2 for line in adata.X.T]

    adata2 = adata[:, adata.var['coverage_feature']>=number_cell_covered].copy()

    for index in range(len(adata2.var_names.tolist())):
        adata2.X[:,index] = np.nan_to_num(adata2.X[:,index], nan=adata2.var['mean_feature_methylation'][index])


    if save!= None:
        adata2.write(save.rstrip('.h5ad')+'.h5ad')
    if copy==False:
        adata = adata2.copy()
    else:
        return(adata2)




def readandimputematrix(file_name, min_coverage=1):
    """
    Temporary function to load and impute methyaltion count matrix into an AnnData object
    
    Parameters
    ----------
    file_name : file name to read and load
    
    min_coverage : minimum number of cells covered for which we keep and impute a variable
    
    Returns
    -------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    
    """
    with open(file_name) as f:
        file = f.readlines()

    # separate annotation from data    
    head_var = file[0]
    head_var = head_var.split('\t')
    # Then, extract the sample names
    sample_names = []
    data_raw = []
    for l in file[1:]:
        l = l.split('\t')
        sample_names.append(l[0])
        data_raw.append(l[1:])

    # clear memory of useless variables 
    del file
    
    ##########################################
    # now, removing empty columns
    empties = []
    partial = []
    full =  []
    for index in range(1, len(data_raw[0])):
        column = [element[index] for element in data_raw]
        if len(list(set(column))) == 1:
            empties.append(index)
        elif len(list(set(column))) <= min_coverage:
            partial.append(index)
        else:
            full.append(index)
         
    ##########################################
    intermed_matrix = []
    name_windows_covered = []
    # let's remove the compltetly uninformative columns
    for index in range(1, len(head_var[1:])):
        if index in full:
            intermed_matrix.append([element[index] for element in data_raw])
            name_windows_covered.append(head_var[index])

    ########################################
    # imputing values.
    imputed_matrix = []
    for row in intermed_matrix:
        imputed_row = []
        if "nan" in row:
            mean = np.mean([float(e)  for e in row if e != "nan"])
            for element in row:
                if element == "nan":
                    imputed_row.append(str(mean))
                else: 
                    imputed_row.append(element)
            imputed_matrix.append(imputed_row)
        else:
            imputed_matrix.append(row)

    imputed_matrix = np.matrix(imputed_matrix).transpose()
    return(ad.AnnData(imputed_matrix, obs=pd.DataFrame(index=sample_names), var=pd.DataFrame(index=name_windows_covered)))
    #return(imputed_matrix, sample_names, name_windows_covered)
