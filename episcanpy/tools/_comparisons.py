import anndata as ad
import numpy as np
import pandas as pd

def transfer_obs(adata1, adata2, obs_key, copy=False):
    """
    adata1 : adata1 receiving the new obs
    adata2 : adata2 providing the additional obs
    obs_key = str or list of str with the obs keys to transfer according to their index
    copy : if True, return a new Anndata otherwise, overwrite adata1
    """
    #extract the data I want to copy
    if isinstance(obs_key, str):
        obs_key = [obs_key]
    markers = adata2.obs[obs_key]
    
    # merging adata1.obs with markers
    mergedDf = adata1.obs.merge(markers,how='outer',
                            left_index=True, right_index=True)

    # then remove the cells that were not in adata4 at the beginning
    annot = []
    for n in mergedDf.index.tolist():
        if n not in adata1.obs_names.tolist():
            annot.append('remove')
        else:
            annot.append('keep')
    mergedDf['filter'] = annot
    mergedDf = mergedDf[mergedDf['filter']=='keep']
    del mergedDf['filter']
    #remove any duplicated lines
    mergedDf = mergedDf[~mergedDf.index.duplicated(keep='first')]
    #reindex:
    mergedDf.index = mergedDf.index.str.strip()
    mergedDf = mergedDf.reindex(adata1.obs_names.tolist())
    mergedDf.index == adata1.obs_names
    
    if copy==True:
        adata3 = adata1.copy()
        adata3.obs = mergedDf.copy()
        return(adata3)
    else:
        adata1.obs = mergedDf.copy()

#redo the imputation for the 3 markers
def remove_values_from_list(the_list, val):
    return([value for value in the_list if value != val])


def imputation(adata, variable_to_input, imputation_variable, var_category=None):
    """
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