import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

def cluster_composition(adata, cluster, condition, xlabel='cell cluster',
                        ylabel='cell count', title=None, save=False):
    """
    """
    

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