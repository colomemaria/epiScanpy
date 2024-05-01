import numpy as np
import pandas as pd

from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from collections import defaultdict

from IPython.display import display
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



def binary_motif_mtx(adata):

    motif_presence_mtx = np.zeros((len(adata.uns["motif_search"]["tf_motifs"]), adata.n_vars))

    for tfbs in adata.uns["motif_search"]["results"]:
        motif_presence_mtx[tfbs["motif_idx"], tfbs["feature_idx"]] = 1

    adata.uns["motif_enrichment"]["motif_presence_mtx"] = motif_presence_mtx



def motif_enrichment(adata, alpha=0.05):

    adata.uns["motif_enrichment"] = {}
    adata.uns["motif_enrichment"]["results"] = defaultdict(list)

    # get motif presence matrix
    binary_motif_mtx(adata)

    results = []

    # get grouping variable from DA test
    grouping_var = adata.uns["rank_features_groups"]["params"]["groupby"]

    M = adata.n_vars  # population size  //  number of peaks
    N = len(adata.uns["rank_features_groups"]["names"])  # sample size  //  number of DA peaks

    # test per motif per group
    for motif_idx in range(len(adata.uns["motif_search"]["tf_motifs"])):

        tf_name = adata.uns["motif_search"]["tf_motifs"][motif_idx].name
        tf_mtx_id = adata.uns["motif_search"]["tf_motifs"][motif_idx].base_id

        n = adata.uns["motif_enrichment"]["motif_presence_mtx"][motif_idx].sum()  # number of successes in the population  //  number of peaks with motif

        for group in adata.obs[grouping_var].cat.categories.tolist():

            da_peaks = adata.uns["rank_features_groups"]["names"][group]
            idcs = [adata.var.index.get_loc(peak) for peak in da_peaks]

            k = adata.uns["motif_enrichment"]["motif_presence_mtx"][motif_idx][idcs].sum()  # number of successes in the sample  // number of DA peaks with motif

            pval = hypergeom.sf(k-1, M, n, N)  # k-1 because hypergeom.sf(k, M, n, N) returns the probability of getting more than k successes
                
            results.append(
                {
                    "tf_name": tf_name,
                    "tf_mtx_id": tf_mtx_id,
                    "pct_tfbs_in_da_features": (k / n) * 100,
                    "pct_da_features_with_tfbs": (k / N) * 100,
                    "pval": pval,
                    "group": group
                }
            )

    pvals = [x["pval"] for x in results]

    _, pvals_corr, _, _ = multipletests(pvals, alpha=alpha, method="fdr_bh")
    for result, pval_corr in zip(results, pvals_corr):
        result["pval_corr_bh"] = pval_corr

    _, pvals_corr, _, _ = multipletests(pvals, alpha=alpha, method="bonferroni")
    for result, pval_corr in zip(results, pvals_corr):
        result["pval_corr_bonferroni"] = pval_corr

    results_grouped = defaultdict(list)
    for result in results:
        results_grouped[result["group"]].append(result)

    # convert defaultdict to dict
    results_grouped = dict(results_grouped)
                
    for group in results_grouped.keys():
        adata.uns["motif_enrichment"]["results"][group] = pd.DataFrame(results_grouped[group])
        adata.uns["motif_enrichment"]["results"][group] = adata.uns["motif_enrichment"]["results"][group].drop(["group"], axis=1)
        adata.uns["motif_enrichment"]["results"][group] = adata.uns["motif_enrichment"]["results"][group].sort_values(by="pval")

    for group in adata.uns["motif_enrichment"]["results"]:
        print(group)
        display(adata.uns["motif_enrichment"]["results"][group])
        print("\n")



def plot_enrichment(x, y, c, s, x_label, colorbar_label, legend_label, show_nonsignificant=True, highlight_significant=True, highlight_list=None, title=None, figsize=(6, 6), dpi=300, save=None):

    if not show_nonsignificant:
        x = x[highlight_list]
        y = y[highlight_list]
        c = c[highlight_list]
        s = s[highlight_list]

    fig, ax = plt.subplots(figsize=figsize, nrows=1, ncols=1, squeeze=True)

    # s_var_scaled = (s / s.max()) * 250
    s_var_scaled = (s - s.min() + 1) / s.max() * 250

    pcol = ax.scatter(
        x=x,
        y=y,
        c=c,
        s=s_var_scaled,
        cmap="viridis_r",
        alpha=[1 if highlight and highlight_significant else 0.4 for highlight in highlight_list] if show_nonsignificant else 1,
        linewidth=1,
        edgecolor="grey",
        vmin=np.floor(c.min()) - 0.25,
        vmax=np.ceil(c.max()) + 0.25
    )

    ax.set_xlabel(x_label, weight="bold")

    plt.yticks(weight="bold")
    
    ax.invert_yaxis()
    ax.yaxis.grid(True, color="lightgrey")
    ax.set_axisbelow(True)

    divider = make_axes_locatable(ax)
    ax_cbar = divider.append_axes(position="right", size=0, pad=1.25)
    ax_cbar.axis("off")
    
    cbar_ticks = np.arange(np.floor(c.min()), np.ceil(c.max()) + 1)

    cbar = plt.colorbar(
        pcol, 
        ax=ax_cbar, 
        shrink=0.625, 
        aspect=10, 
        pad=0, 
        anchor=(0, 0.05),
        ticks=cbar_ticks
    )
    
    cbar.ax.set_yticklabels(["$10^{{{:n}}}$".format(x) for x in cbar_ticks])
    cbar.set_label(label=colorbar_label, weight="bold", labelpad=12)

    legend_size_vals_min = s.min()
    legend_size_vals_max = s.max()
    legend_size_vals = np.linspace(start=legend_size_vals_min, stop=legend_size_vals_max, num=4)
    legend_size_vals = np.round(legend_size_vals)
    # legend_size_vals_scaled = (legend_size_vals / legend_size_vals_max) * 250
    legend_size_vals_scaled = (legend_size_vals - legend_size_vals.min() + 1) / legend_size_vals_max * 250

    legend = ax.legend(
        [plt.scatter([], [], s=i, color="black") for i in legend_size_vals_scaled],
        ["{:n}".format(i) for i in legend_size_vals],
        title=legend_label,
        title_fontproperties={"weight": "bold", "size": "small"},
        fontsize="small",
        loc="upper left",
        bbox_to_anchor=(1.05 ,1),
        frameon=True,
        labelspacing=1,
        borderpad=0.8
    )

    ylim = ax.get_ylim()
    ylim = (ylim[0] - (ylim[1] - ylim[0]) * 0.05, ylim[1] + (ylim[1] - ylim[0]) * 0.05)
    ax.set_ylim(ylim)

    legend.get_title().set_horizontalalignment("center")

    if title is not None:
        fig.suptitle(title, weight="bold", y=0.95)
    
    if save is None:
        plt.show()
    else:
        plt.savefig(save, dpi=dpi, bbox_inches="tight")



def plot_motif_enrichment(adata, corr_method="benjamini-hochberg", top_n=10, show_nonsignificant=True, highlight_significant=True, pval_threshold=0.05, dpi=300, save=None):

    if corr_method is None or corr_method.lower() in ["raw", "none"]:
        pval_key = "pval"
        colorbar_label = "unadj. p-value"
    
    elif corr_method.lower() in ["benjamini-hochberg", "bh"]:
        pval_key = "pval_corr_bh"
        colorbar_label = "adj. p-value (Benjamini-Hochberg)"
    
    elif corr_method.lower() in ["bonferroni"]:
        pval_key = "pval_corr_bonferroni"
        colorbar_label = "adj. p-value (Bonferroni)"

    msgs = []

    for i, group in enumerate(adata.uns["motif_enrichment"]["results"]):

        data = adata.uns["motif_enrichment"]["results"][group]

        x = data.pct_tfbs_in_da_features
        y = data.tf_name + " (" + data.tf_mtx_id + ")"
        c = data[pval_key]
        s = data.pct_da_features_with_tfbs

        if top_n:
            x, y, c, s = x[:top_n], y[:top_n], c[:top_n], s[:top_n]

        highlight_list = [True if pval < pval_threshold else False for pval in c]

        if not show_nonsignificant and  np.sum(highlight_list) == 0:
            msgs.append(f"No significant results for {group}\n")
            continue

        plot_enrichment(
            x, 
            y, 
            np.log10(c), 
            s, 
            x_label="% TFBS in DA features", 
            colorbar_label=colorbar_label, 
            legend_label="% DA features\nwith TFBS", 
            show_nonsignificant=show_nonsignificant,
            highlight_significant=highlight_significant,
            highlight_list=highlight_list,
            title=group, 
            figsize=(8, 6),
            dpi=dpi,
            save=save.split(".")[0] + f"_{i+1}." + save.split(".")[1] if save is not None else None
        )

    for msg in msgs:
        print(msg)