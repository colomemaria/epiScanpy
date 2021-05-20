import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad
import pandas as pd

import warnings
from warnings import warn

from scipy.sparse import issparse
from scipy.stats.stats import pearsonr, spearmanr

def _correlation(adata, variable, component, use_rep, method):
    x = adata.obsm[use_rep][:,component-1]
    y = adata.obs[variable]
        
    if method=='pearson':
        correlation = pearsonr(x, y)
    elif method=='spearman':
        correlation = spearmanr(x, y)
        
    return(correlation)

def correlation_component(adata,
                   variable,
                   component='all',
                   use_rep='X_pca',
                   method='pearson',
                   absolute = True,
                   title=None,
                   xlabel=None,
                   ylabel=None,
                   color_palette=None,
                   #figsize=(15,10),
                   #save_dpi=250,
                   show=True,
                   save=None,
                  ):
    """
    Correlation between a given PC and a covariate.     
    If show == True, plot a scatter plot. 
    Available methods for correlation: 'pearson' or 'spearman'
    
    Parameters
    ----------
    
    adata : input adata
    
    variable : covariate either saved in obs or var
    
    component : 'all', int or list of int corresponding to the components to use. 
    start at 1st component. Do not specify 0 as 1st component.
    if all, compute correlation for all components
    
    use_rep : 'X_pca', 'X_fa', 'X_nmf', 'X_lsi'. depends of the decomposition used. 
    
    method : 'pearson' or 'spearman' available - scipy.stats implementation
    
    absolute : bool if True, return absolute correlation value 
    
    color_palette : seaborn color palette to use for plotting
    
    title : optional title to the plot
    
    xlabel : optional xlabel
    
    ylabel : optional ylabel
    
    show: Print the correlation coefficient and p-value. Additionaly, plot a scatter plot
    of the PC coordinate and the covariate value for every cell or feature in your matrix
    
    save: if specified the name of the picture, save as png
    
    Return 
    ------
    
    correlation coefficient and p-value
    
    """
    all_correlation={'correlation':[], 'pval':[]}
    if isinstance(component, int):
        component = [component]
    elif component == 'all':
        component = list(range(1, adata.obsm[use_rep].shape[1]+1,1))
        
    if np.max(component)-1>len(adata.varm['PCs'][0]):
        warnings.warn("".join(["""You requested a component that is not currently available.
                            If you used X_pca decomposition, please run epi.pp.pca(adata, n_comps=""", int(np.max(component)+1), ') ']))
        
    for n in component:
        correlation = _correlation(adata=adata,
                                   variable=variable,
                                   component=n,
                                   use_rep=use_rep,
                                   method=method)
        all_correlation['correlation'].append(correlation[0])
        all_correlation['pval'].append(correlation[1])
    
    tmp_list = []
    for value in all_correlation['correlation']:
        if value < 0:
            tmp_list.append(-1*value)
        else:
            tmp_list.append(value)
    all_correlation['abs_correlation'] = tmp_list
    
    adata.uns['correlation_components'] = all_correlation
    
    #fig = plt.figure(figsize=figsize)
    if show and len(component)==1:
        component = component[0]
        x = adata.obsm[use_rep][:,component-1]
        y = adata.obs[variable]
        ## Add legends and title
        if xlabel:
            plt.xlabel(xlabel)
        else:
            plt.xlabel(" ".join([use_rep, str(component)]))
            
        if ylabel:
            plt.ylabel(ylabel)
        else:
            plt.ylabel(variable)
            
        if title:
            plt.title(title)
            
        plt.scatter(x, y)
        print("correlation: " +str(correlation[0]))
        print("pval: " +str(correlation[1]))
        
    elif show:
        _plot_correlations(all_correlation,
                           absolute=absolute,
                           color_palette=color_palette,
                           save=None)
        
    if save!= None:
        #fig.savefig(save, dpi=save_dpi)
        plt.savefig(save, bbox_inches="tight")
    plt.show()
    
def _plot_correlations(all_correlation, absolute=True, color_palette=None, save=None):
    df = pd.DataFrame.from_dict(all_correlation)
    df['components'] = [x+1 for x in df.index]
    if color_palette == None:
        if absolute == True:
            color_palette="viridis"
        else:
            color_palette='icefire'
    
    if absolute == True:
        del df['pval'], df['correlation']
        g = sns.PairGrid(df.sort_values("components", ascending=False),
                     x_vars=df.columns[:-1], y_vars=['components'],
                     height=5+len(df['abs_correlation'])/5, palette=color_palette)
        # Draw a dot plot using the stripplot function
        g.map(sns.stripplot, size=10, orient="h", jitter=False,
              linewidth=1, edgecolor="w")

        # Use the same x axis limits on all columns and add better labels
        g.set(xlim=(-0.1, 1), xlabel="Absolute correlation value per component", ylabel="")
        
    else:
        del df['pval'], df['abs_correlation']
        g = sns.PairGrid(df.sort_values("components", ascending=False),
                     x_vars=df.columns[:-1], y_vars=['components'],
                     height=5 +len(df['correlation'])/5, palette=color_palette)
        # Draw a dot plot using the stripplot function
        g.map(sns.stripplot, size=10, orient="h", jitter=False,
              linewidth=1, edgecolor="w")

        # Use the same x axis limits on all columns and add better labels
        g.set(xlim=(-1, 1), xlabel="Correlation value per component", ylabel="")

    if save!=None:
        plt.savefig(save)
    plt.show()

def filtering_components(adata, components, use_rep='X_pca', new_rep=None):
    """
    Filter out components.
    
    Parameters
    ----------
    adata : AnnData object
    
    components : int or list of int corresponding to the PCs to remove. 
    
    use_rep : 'X_pca', 'X_lsi' etc. obsm key
    
    new_rep : if None, overwrite use_rep key. Else, use new_rep as a new obsm key. 
    """
    if isinstance(components, int):
        components=[components]
    if new_rep == None:
        new_rep = use_rep
    
    components = [x-1 for x in components]
    adata.obsm[new_rep] =  np.delete(adata.obsm[use_rep], components, axis=1)