import numpy as np
from scipy.stats import hypergeom


def binary_motif_mtx(adata):

    motif_presence_mtx = np.zeros((len(adata.uns["motif_search"]["tf_motifs"]), adata.n_vars))

    for tfbs in adata.uns["motif_search"]["results"]:
        motif_presence_mtx[tfbs["motif_idx"], tfbs["feature_idx"]] = 1

    adata.uns["motif_search"]["motif_presence_mtx"] = motif_presence_mtx