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

def imputation_met(adata, number_cell_covered=10, imputation_value='mean', save=None):
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

    # This step need to be sped up and coule be multithread.
    # Only the mean available for now. And only the minimum number of cells covered and not the variety of the 
    # methylation levels
    # also, it odes not return the variable annoations and force to add 2 values
    matrix = adata.X.transpose().tolist()
    old_features = adata.var_names.tolist()
    
    new_matrix = []
    new_features_name = []
    means = []
    medians = []
    feat_nb = 0
    for feature in adata.X:
        nb_cov = 0
        # also the different number of values
        val_for_mean = []
        for f in feature:
            if 0.0 <= f <= 1.0:
                nb_cov += 1
                val_for_mean.append(f)
    
    
        new_feature = []
    
        if nb_cov >= number_cell_covered:
            new_features_name.append(old_features[feat_nb])
            # also median possible
            mean = np.mean(val_for_mean)
            for f in feature:
                if 0.0 <= f <= 1.0:
                    new_feature.append(f)
                else:
                    new_feature.append(mean)
                
            new_matrix.append(new_feature)
        
            means.append(mean)
            medians.append(np.median(val_for_mean))
            
        feat_nb += 1
    
    
    new_matrix = np.array(new_matrix)
    return(new_matrix, new_features_name)
    adata2 = ad.AnnData(new_matrix.transpose(), obs=adata.obs_names, var = pd.DataFrame(index=new_features_name))
    adata2.obs = adata.obs
    adata2.var['median_before_imput'] = medians # before imputation
    adata2.var['mean_before_imput'] = means # before imputation
    
    if save == None:
        return(adata2)
    else:
        if (len(save)>=5) and (save[-5:] == '.h5ad'):
            adata2.write(save)
        else:
            adata2.write(save+'.h5ad')
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
    return(imputed_matrix, sample_names, name_windows_covered)
