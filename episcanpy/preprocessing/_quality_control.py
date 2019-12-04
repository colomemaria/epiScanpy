import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

import warnings
from warnings import warn

from scipy.sparse import issparse
from scipy.stats.stats import pearsonr, spearmanr

def cal_var(adata, show=True):
    
    if issparse(adata.X):
        adata.var['n_cells'] = adata.X.sum(axis=0).tolist()[0] ## double check this
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = abs(adata.var['prop_shared_cells']-0.5)
    else:
        adata.var['n_cells'] = adata.X.sum(axis=0)
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = abs(adata.var['prop_shared_cells']-0.5)
    
    if show: # plotting
        
        fig = plt.figure(figsize=(8,6))
        sns.distplot(adata.var['variability_score'], bins=40)
        sns.distplot(adata.var['prop_shared_cells'], bins=40)
        fig.legend(labels=['variability_score','prop_shared_cells'],
                   loc='upper right')
        plt.ylabel('nb of features')
        plt.title('Distribution of feature coverage')
        plt.show()
        
def select_var_feature(adata, max_score=0.5, nb_features=None, show=True, copy=False):
    """
    This function computes a variability score to rank the most variable features across all cells. 
    Then it selects the most variable features according to either a specified number of features (nb_features) or a maximum variance score (max_score).
    
    Parameters
    ----------

    adata: adata object
    
    max_score: max threshold of the variability score to retain features,
    where 0 is the score of the most variable features and 0.5 is the score of the least variable features.
    
    nb_features: default value is None, if specify it will select a the top most variable features.
    if the nb_features is larger than the total number of feature, it filters based on the max_score argument
    
    show: default value True, it will plot the distribution of var.
    
    copy: return a new adata object if copy == True.

    Returns
    -------
    Depending on ``copy``, returns a new AnnData object or overwrite the input


    """
    if copy:
        inplace=False
    else:
        inplace=True

    adata = adata.copy() if not inplace else adata
    
    # calculate variability score
    cal_var(adata, show=show) # adds variability score for each feature 
    # adata.var['variablility_score'] = abs(adata.var['prop_shared_cells']-0.5)
    var_annot = adata.var.sort_values(ascending=True, by ='variability_score')

    # calculate the max score to get a specific number of feature        
    if nb_features != None and nb_features < len(adata.var_names):  
            max_score = var_annot['variability_score'][nb_features]    
     
    
    adata_tmp = adata[:,adata.var['variability_score']<=max_score].copy()
    
    ## return the filtered AnnData objet.
    if not inplace:
        adata_tmp = adata[:,adata.var['variability_score']<=max_score]
        return(adata_tmp)
    else:
        adata._inplace_subset_var(adata.var['variability_score']<=max_score)
        
def binarize(adata, copy=False):
    """convert the count matrix into a binary matrix.

    Parameters
    ----------

    adata: AnnData object

    copy: return a new adata object if copy == True.

    Returns
    -------
    Depending on ``copy``, returns a new AnnData object or overwrite the input

    """

    if copy:
        adata2 = adata.copy()
        adata2.X[adata2.X != 0] = 1
        return(adata2)
    else:
        adata.X[adata.X != 0] = 1
        

def commonness_features(adata,
                        binary=None, 
                        log=False,
                        key_added=None, 
                        
                        threshold=None,
                        bw=0.5,
                        bins=50,
                        xlabel=None, ylabel=None, title=None,
                        color=None, edgecolor=None,
                        
                        save=None):
    
    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.

    Parameters
    ----------
    adata
    threshold
    log
    binary
    key_added
    bw
    xlabel
    title
    color
    edgecolor
    
    """
    
    if key_added == None:
        key_added='commonness'
        
    #if xlabel ==None:
    #    plt.xlabel('cells sharing a feature')
    #else:
    #    plt.xlabel(xlabel)    

    # plotting parameters
    if xlabel ==None:
        plt.xlabel('number of features')
    else:
        plt.xlabel(xlabel)
        
    if ylabel ==None:
        plt.ylabel('number of cells')
    else:
        plt.ylabel(ylabel)
    
    if title !=None:
        plt.title(title)
    
    
    # calculate for each feature in how many cell it is open
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        common = np.sum(adata.X, axis=0).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        common = np.sum(tmp_array.X, axis=0).tolist()

    if issparse(adata.X):
        common = common[0]        

    # save key
    adata.var[key_added] = common

    if threshold != None:
        plt.axvline(x=threshold, color='r')
    bw_param = bw
    sns.set_style('whitegrid')
    #sns.kdeplot(np.array(adata.var[key_added]), bw=bw_param)


    if color == None:   
        color='c'
    if edgecolor == None:
        edgecolor='k'

    if log !=False:
        if 0 in common:
            warnings.warn("""Some features are not open in any cell. Use epi.pp.filter_features(adata, max_cells=1) to remove these features.""")
        
        # different potential bases for log transfo
        if log=='log2':
            plt.xlabel('number of features (log2)')
            fig = plt.hist(np.log2(common), bins, color=color, edgecolor=edgecolor)
        elif log=='log1p':
            plt.xlabel('number of features (log1p)')
            fig = plt.hist(np.log1p(common), bins, color=color, edgecolor=edgecolor)
        elif log=='log':
            plt.xlabel('number of features (log natural base e)')
            fig = plt.hist(np.log(common), bins, color=color, edgecolor=edgecolor)
        else:
            plt.xlabel('number of features (log10)')
            fig = plt.hist(np.log10(common), bins, color=color, edgecolor=edgecolor)
    
    else:
        fig = plt.hist(common, bins=int(80), color=color, edgecolor=edgecolor)
 
    
    if title !=None:
        plt.title(title)
    
    #fig = plot.get_figure()
    if save!= None:
        fig.savefig(save)
    plt.show()

    adata.var[key_added] = common


def coverage_cells(adata,
                   key_added=None,
                   log=False,
                   binary=None,
                   ## threshold 
                   threshold=None,
                   bw=0.5,
                   ## plotting
                    bins=50,
                   xlabel=None, ylabel=None, title=None,
                   color=None, edgecolor=None,
                   ## sving
                   save=None):
    """
    Histogram of the number of open features (in the case of ATAC-seq data) per cell.

    Parameters
    ----------

    adata
    
    binary
    log : default not . Accepted log2, log10, log1p, log (natural log base e).
    if log is True, log10 used

    bins
    threshold

    key_added

    xlabel
    ylabel
    title
    color
    edgecolor
    save

    """
    if key_added == None:
        key_added='nb_features'
    
    # calculate the number of features per cell
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        sum_peaks = np.sum(adata.X, axis=1).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        sum_peaks = np.sum(tmp_array.X, axis=1).tolist()
  
    
    if issparse(adata.X):
        sum_peaks = [element[0] for element in sum_peaks]        

    adata.obs[key_added] = sum_peaks
    
    
    # plotting parameters
    if xlabel ==None:
        plt.xlabel('number of features')
    else:
        plt.xlabel(xlabel)
        
    if ylabel ==None:
        plt.ylabel('number of cells')
    else:
        plt.ylabel(ylabel)
    
    if title !=None:
        plt.title(title)
        
    if color == None:   
        color='c'
    if edgecolor == None:
        edgecolor='k'

    if threshold != None:
        plt.axvline(x=threshold, color='r')
    bw_param = bw

    if log!=False:
        if 0 in sum_peaks:
            warnings.warn("""Some cells do not contain any open feature. Use epi.pp.filter_cells(adata, max_features=1) to remove these cells.""")
        
        if log=='log2':
            plt.xlabel('number of features (log2)')
            fig = plt.hist(np.log2(sum_peaks), bins, color=color, edgecolor=edgecolor)
        elif log=='log1p':
            plt.xlabel('number of features (log1p)')
            fig = plt.hist(np.log1p(sum_peaks), bins, color=color, edgecolor=edgecolor)
        elif log=='log':
            plt.xlabel('number of features (log natural base e)')
            fig = plt.hist(np.log(sum_peaks), bins, color=color, edgecolor=edgecolor)
        else:
            plt.xlabel('number of features (log10)')
            fig = plt.hist(np.log10(sum_peaks), bins, color=color, edgecolor=edgecolor)

    else:
        fig = plt.hist(sum_peaks, bins, color=color, edgecolor=edgecolor)
    
    #fig = plot.get_figure()
    if save!= None:
        plt.savefig(save)

    plt.show()
    adata.obs[key_added] = sum_peaks


def correlation_PC(adata,
                   variable,
                   pc=1,
                   obs=True,
                   method='pearson',
                   title=None.
                   xlabel=None,
                   ylabel=None,
                   show=True,
                   save=None,
                  ):
    """
    Correlation between a given PC and a covariate.     
    If show == True, plot a scatter plot. 
    If obs == True, consider a obs covariate. Else, take a variable covariate.
    Available methods for correlation: 'pearson' or 'spearman'
    
    Parameters
    ----------
    
    adata : input adata
    
    variable : covariate either saved in obs or var
    
    obs : Boolean that specify if the covariate is stored in obs or in var.
    In later version, this parameter will not need to be specified
    
    method : 'pearson' or 'spearman' available - scipy.stats implementation
    
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
    if pc-1>len(adata.varm['PCs'][0]):
        warnings.warn("".join(["""You requested a PC that is not currently available.
                            Please run epi.pp.pca(adata, n_comps=""", int(pc+1), ') ']))
        
    
    if obs:
        x = adatareduced.obsm['X_pca'][:,pc-1]
        y = adatareduced.obs[variable]
    else:
        x = adatareduced.varm['PCs'][:,pc-1]
        y = adatareduced.var[variable]
        
    if method=='pearson':
        correlation = pearsonr(x, y)
    elif method=='spearman':
        correlation = spearmanr(x, y)
    
    if show:
        ## Add legends and title
        if xlabel:
            plt.xlabel(xlabel)
        else:
            plt.xlabel(" ".join(['X_PC', pc]))
            
        if ylabel:
            plt.ylabel(ylabel)
        else:
            plt.ylabel(variable)
            
        if title:
            plt.title(title)
            
        plt.scatter(x, y)
        print(correlation)
        
        
    if save!= None:
        plt.savefig(save)
    
    
    return (correlation)
    

