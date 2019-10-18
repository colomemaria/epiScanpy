import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax

# function to calculate the variance
# function to calculate the variance
def cal_var(adata, show=True):
    
    if str(type(adata.X)) != "numpy.ndarray":
        adata.var['n_cells'] = adata.X.sum(axis=0)
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variablility_score'] = abs(adata.var['prop_shared_cells']-0.5)
    else:
        adata.var['n_cells'] = adata.X.sum(axis=0).tolist()[0]
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variablility_score'] = abs(adata.var['prop_shared_cells']-0.5)
    
    if show: # plotting
        
        fig = plt.figure(figsize=(8,6))
        sns.distplot(adata.var['variablility_score'], bins=40)
        sns.distplot(adata.var['prop_shared_cells'], bins=40)
        fig.legend(labels=['variablility_score','prop_shared_cells'],
                   loc='upper right')
        plt.ylabel('nb of features')
        plt.title('Distribution of feature coverage')
        plt.show()