import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from collections import defaultdict


def binary_motif_mtx(adata):

    motif_presence_mtx = np.zeros((len(adata.uns["motif_search"]["tf_motifs"]), adata.n_vars))

    for tfbs in adata.uns["motif_search"]["results"]:
        motif_presence_mtx[tfbs["motif_idx"], tfbs["feature_idx"]] = 1

    adata.uns["motif_enrichment"]["motif_presence_mtx"] = motif_presence_mtx


def motif_enrichment(adata, pval_threshold=0.05):

    adata.uns["motif_enrichment"] = {}
    adata.uns["motif_enrichment"]["results"] = defaultdict(list)

    # get motif presence matrix
    binary_motif_mtx(adata)

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

            if pval < pval_threshold:
                
                adata.uns["motif_enrichment"]["results"][group].append(
                    {
                        "tf_name": tf_name,
                        "tf_mtx_id": tf_mtx_id,
                        "pval": pval
                    }
                )

    # convert defaultdict to dict
    adata.uns["motif_enrichment"]["results"] = dict(adata.uns["motif_enrichment"]["results"])
                
    for group in adata.uns["motif_enrichment"]["results"]:
        adata.uns["motif_enrichment"]["results"][group] = pd.DataFrame(adata.uns["motif_enrichment"]["results"][group])

    # sort results by p-value
    for group in adata.uns["motif_enrichment"]["results"]:
        adata.uns["motif_enrichment"]["results"][group] = adata.uns["motif_enrichment"]["results"][group].sort_values(by="pval")

    # sort results by p-value
    for group in adata.uns["motif_enrichment"]["results"]:
        print(group)
        print(adata.uns["motif_enrichment"]["results"][group], end="\n\n")
