import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
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

    # create dataframe
    df = pd.crosstab(adata.obs[obs_1], adata.obs[obs_2])
    array = np.array(df)
    x = df.columns.tolist()
    y = df.index.tolist()

    #colors for the plot
    if obs_1+"_colors" in adata.uns.keys():
        colors=adata.uns[obs_1+"_colors"]
    else:
        # select random colors
        no_of_colors=len(y)
        colors=["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
           for j in range(no_of_colors)]

    # plot bars in stack manner
    previous_value = 0
    index = 0
    for n in range(len(y)):
        plt.bar(x, array[index], bottom=previous_value, color=colors[index])
        previous_value += array[index]
        index += 1

    #    The strings
    #    ``'upper left', 'upper right', 'lower left', 'lower right'``
    #    place the legend at the corresponding corner of the axes/figure.
    #    The strings
    #    ``'upper center', 'lower center', 'center left', 'center right'``
    #    place the legend at the center of the corresponding edge of the
    #    axes/figure.
    plt.xticks(x, rotation=90)
    plt.legend(y, loc=loc_legend, bbox_to_anchor=location_bbox)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if save != None :
        plt.savefig(save, bbox_inches='tight')
    plt.show()


def histogram(adata,
              key,
              bins=40,
              min_threshold=None,
              max_threshold=None,
              show_log=True,
              show_mean=True,
              show_median=True,
              print_statistics=True,
              figsize=None,
              save=None):
    """
    General purpose histogram plotting function.

    Args:
        adata: AnnData object
        key: name of the variable. has to be either in .obs or .var
        bins: number of bins
        min_threshold: minimum threshold
        max_threshold: maximum threshold
        show_log: show logarithmized (log10) data
        show_mean: show mean (dashed line)
        show_median: show median (straight line)
        print_statistics: print basic statistics
        figsize: size of the figure
        save: if True or str, save the figure. str represents entire path. filetype is inferred.

    Returns:
        None
    """

    if figsize is None:
        figsize = (6, 4) if not show_log else (12, 4)

    ncols = 1 if not show_log else 2

    fig, axs = plt.subplots(figsize=figsize, nrows=1, ncols=ncols, squeeze=False)
    axs = axs.flatten()

    df = adata.obs if key in adata.obs else adata.var

    axs[0].hist(df[key], bins=bins, linewidth=0.5, edgecolor="black")

    if min_threshold:
        axs[0].axvline(x=min_threshold, color="red", linestyle="--", linewidth=1, alpha=0.75)
    if max_threshold:
        axs[0].axvline(x=max_threshold, color="red", linestyle="--", linewidth=1, alpha=0.75)

    if show_mean:
        axs[0].axvline(x=df[key].mean(), color="white", linestyle="--", linewidth=2, alpha=0.75)
    if show_median:
        axs[0].axvline(x=df[key].median(), color="white", linestyle="-", linewidth=2, alpha=0.75)

    axs[0].set_xlabel("{}".format(key))

    if show_log:
        axs[1].hist([np.log10(val) if val != 0 else 0 for val in df[key].values], bins=bins, linewidth=0.5, edgecolor="black")

        if min_threshold:
            axs[1].axvline(x=np.log10(min_threshold), color="red", linestyle="--", linewidth=1, alpha=0.75)
        if max_threshold:
            axs[1].axvline(x=np.log10(max_threshold), color="red", linestyle="--", linewidth=1, alpha=0.75)

        if show_mean:
            axs[1].axvline(x=np.log10(df[key].mean()), color="white", linestyle="--", linewidth=2, alpha=0.75)
        if show_median:
            axs[1].axvline(x=np.log10(df[key].median()), color="white", linestyle="-", linewidth=2, alpha=0.75)

        axs[1].set_xlabel("log10({})".format(key))

    plt.tight_layout()

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "qc_histogram_{}.png".format(key)

        plt.savefig(filename, dpi=300)

    if print_statistics:
        print("Max:\t{}".format(df[key].max()))
        print("Median:\t{}".format(df[key].median()))
        print("Mean:\t{}".format(df[key].mean()))
        print("Min:\t{}".format(df[key].min()))
