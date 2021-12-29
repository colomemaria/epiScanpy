import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

# importing package
import matplotlib.pyplot as plt
import numpy as np
import random

import warnings
from warnings import warn

def cluster_composition(adata, cluster, condition, xlabel='cell cluster',
                        ylabel='cell count', title=None, save=False):
    """
    Deprecated. Use epi.pl.cell_composition instead.
    """
    warnings.warn("Deprecated. Use epi.pl.cell_composition instead.")

    contingency_table = pd.crosstab(
        adata.obs[condition],
        adata.obs[cluster],
        margins = True
    )

    counts = []
    p_part = []
    index = 0
    categories = sorted(list(set(adata.obs[cluster])))
    for n in sorted(set(adata.obs[condition])):
        #counts.append()
        p_part.append(plt.bar(categories, contingency_table.iloc[index][0:-1].values))
        index += 1

    #Plots the bar chart
    #plt.figsize(figsize=[6.4, 4.8])
    plt.legend(tuple([p[0] for p in p_part]), tuple(sorted(set(adata.obs[condition]))))
    plt.xlabel(xlabel, )
    plt.ylabel(ylabel)
    plt.title(title)
    
    
    if save!=False:
        
        if (save==True) or (save.split('.')[-1] not in ['png', 'pdf']):
            plt.savefig('cluster_composition.png', dpi=300, bbox_inches="tight")
        else:
            plt.savefig('_'.join(['cluster_composition',save]), #format=save.split('.')[-1],
                        dpi=300, bbox_inches="tight")
            
    plt.show()


def cell_composition(adata, obs_1,  obs_2,
                    title='Cell composition per sample',
                    xlabel="",
                    ylabel="Number of cells",
                    loc_legend = 'best',
                    location_bbox=(1, 0, 0, 1),
                    save=None):
    
    """
    Bar plots displaying the cell composition division between two Anndata obs categories. 
    
    adata : AnnData objct
    obs_1 : adata.obs key 1 
    obs_2 : adata.obs key 2
    title : [optional] title of the plot
    loc_legend : location of the legend. Available are ``'upper left', 'upper right', 'lower left', 'lower right'``
    or ``'upper center', 'lower center', 'center left', 'center right'``
    bbox_to_anchor : tuple containing the location of the figure. Default (1, 0, 0, 1)
    save : if not None, str corresponding to the output file name
    """
    
    
    #colors for the plot
    if obs_2+"_colors" in adata.uns.keys():
        colors=adata.uns[obs_2+"_colors"]
    else:
        # select random colors
        no_of_colors=len(df.index.tolist())
        colors=["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
           for j in range(no_of_colors)]


    # create dataframe
    df = pd.crosstab(adata.obs[obs_1], adata.obs[obs_2])
    array = np.array(df)
  
    # plot bars in stack manner
    previous_value = 0
    index = 0
    for n in range(len(df.index.tolist())):
        plt.bar(x, array[index], bottom=previous_value, color=colors[index])
        previous_value +=array[index]
        index += 1

    #    The strings
    #    ``'upper left', 'upper right', 'lower left', 'lower right'``
    #    place the legend at the corresponding corner of the axes/figure.
    #    The strings
    #    ``'upper center', 'lower center', 'center left', 'center right'``
    #    place the legend at the center of the corresponding edge of the
    #    axes/figure.
    plt.xticks(x, rotation=90)
    plt.legend(df.index.tolist(), loc=loc_legend, bbox_to_anchor=location_bbox)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if save != None :
        plt.savefig(save, bbox_inches='tight')
    plt.show()