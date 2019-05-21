import matplotlib.pyplot as plt
import anndata as ad 
import scanpy.api as sc
from sklearn.metrics import silhouette_score, silhouette_samples
import matplotlib.cm as cm
import numpy as np

def silhouette(adata_name, cluster_annot, key=None,
              xlabel=None, ylabel=None, title=None, size='large',
               name_cluster=True, name_cluster_pos='left', 
              palette=None, save=None):
    """
    Plot the product of tl.silhouette as a silhouette plot

    Parameters
    ----------

    adata_name: AnnData object

    cluster_annot: observational variable corresponding to a cell clustering

    key: specify name of precomputed silhouette scores if not standard


    Return
    ------

    Silhouette plot

    """

    
    silhouette_avg = silhouette_score(X, cluster_labels, metric)
    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric)
    
    if key!=None:
        sample_silhouette_values = adata_name.obs[key]
        silhouette_avg = adata_name.uns[key]
    else:
        adata_name.obs['silhouette_samples'] = sample_silhouette_values
        sample_silhouette_values = adata_name.obs['silhouette_samples']
        silhouette_avg = adata_name.obs['silhouette_samples']
        
    
    fig, ax1 = plt.subplots(1, 1)
    if size=='large':
        fig.set_size_inches(18, 7)
    #ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    y_lower = 10
    
    ## for title argument
    if title != None:
        ax1.set_title(title)
    else:
        ax1.set_title(cluster_annot)
    ## axis label arguments
    if xlabel != None:
        ax1.set_xlabel(xlabel)
    if ylabel == None:
        ax1.set_ylabel(cluster_annot)
    else:
        ax1.set_ylabel(ylabel)

        
    i =0
    for name in list(sorted(set(adata_name.obs[cluster_annot]))):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == name]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        
        if palette != None:
            color = palette[i]
        else:
            color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                            0, ith_cluster_silhouette_values,
                            facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        #ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        if name_cluster:
            if name_cluster_pos == 'left':
                ax1.text((min(sample_silhouette_values)+0.1), y_lower + 0.5 * size_cluster_i, name)
            else:
                ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, name)
            

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
        i += 1


    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    #ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    if save != None:
        plt.savefig('_'.join(['silhouette', save])) 
    plt.show()
    
    print(silhouette_avg)
    return()

def silhouette_tot(adata_name, cluster_annot, value='X_pca', metric='euclidean',
              xlabel=None, ylabel=None, title=None, size='large',
               name_cluster=True, name_cluster_pos='left', 
              palette=None, save=None, key_added=None):
    """
    Both compute silhouette scores and plot it.

    It computes the general silhouette score as well as a silhouette score for every cell according 
    to the cell cluster assigned to it. 

    Parameters
    ----------
    adata_name: AnnData object

    cluster_annot: observational variable corresponding to a cell clustering

    value: measure used to build the silhouette plot (X_pca, X_tsne, X_umap)

    metric: 'euclidean'

    key_added: key to save the computed silhouette scores

    Return
    ------

    general silhouette score in 'uns' of the AnnData object
    individual silhouette scores in 'obs' of the AnnData object
    
    Silhouette plot


    Credit to sklearn script : 
    https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html#sphx-glr-auto-examples-cluster-plot-kmeans-silhouette-analysis-py
    return score and silhouette plot. Still some work to do to finish the function.
    size=None but you can put 'large' if you want a bigger default figure size
    """
    #sc.pl.umap(adata_name, color=cluster_annot)
    #print(list(sorted(set(adata_name.obs[cluster_annot]))))

    X = adata_name.obsm[value]
    cluster_labels = adata_name.obs[cluster_annot]
    n_clusters = len(set(adata_name.obs[cluster_annot]))


    fig, ax1 = plt.subplots(1, 1)
    if size=='large':
        fig.set_size_inches(18, 7)
    #ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    y_lower = 10
    
    ## for title argument
    if title != None:
        ax1.set_title(title)
    else:
        ax1.set_title(cluster_annot)
    ## axis label arguments
    if xlabel != None:
        ax1.set_xlabel(xlabel)
    if ylabel == None:
        ax1.set_ylabel(cluster_annot)
    else:
        ax1.set_ylabel(ylabel)

    ## depending of verbosity, print(silhouette_avg)
    ## also, return sample_silhouette_values as adata.obs['silhouette_samples']
    silhouette_avg = silhouette_score(X, cluster_labels, metric)
    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric)
    if key_added:
        adata_name.obs[key_added] = sample_silhouette_values
    else:
        adata_name.obs['silhouette_samples'] = sample_silhouette_values


    i =0
    for name in list(sorted(set(adata_name.obs[cluster_annot]))):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == name]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        
        if palette != None:
            color = palette[i]
        else:
            color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                            0, ith_cluster_silhouette_values,
                            facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        #ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        if name_cluster:
            if name_cluster_pos == 'left':
                ax1.text((min(sample_silhouette_values)+0.1), y_lower + 0.5 * size_cluster_i, name)
            else:
                ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, name)
            

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
        i += 1


    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    #ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    if save != None:
        plt.savefig('_'.join(['silhouette', save])) 
    plt.show()
    
    print(silhouette_avg)
    return()

