import numpy as np
import anndata as ad

def remove_values_from_list(the_list, val):
    return([value for value in the_list if value != val])


def imputation_feature(adata, variable_to_input, imputation_variable, var_category=None):

    """
    impute the missing methylation feature according to the clusters of interest.
    """

    if var_category:
        col_index = adata.var[var_category].index(variable_to_input)
    else:
        col_index = adata.var_names.tolist().index(variable_to_input)
    X = adata.X[:, [col_index]]
    guess_imputed_value = np.median(X)


    # get the new imputer value per louvain cluster:
    new_imputed = {}
    for cluster in list(set(adata.obs[imputation_variable])):
        tmp_adata = adata[adata.obs[imputation_variable]==cluster,:].copy()
        tmp_X = tmp_adata.X[:, [col_index]]
        tmp_X = remove_values_from_list(tmp_X, guess_imputed_value)
        new_imputed[cluster] = np.mean(tmp_X)
    
    annot = []
    annot2 = []
    index = 0
    for cluster in adata.obs[imputation_variable]:
        if X[index] == guess_imputed_value:
            annot.append(new_imputed[cluster])
            annot2.append(np.nan)
        else:
            annot.append(X[index][0])
            annot2.append(X[index][0])
        index +=1 
    adata.obs[variable_to_input+'_no_input'] = annot2
    adata.obs[variable_to_input+'_imputed'] = annot