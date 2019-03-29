import scanpy.api as sc

def plot_rank_features(adata, groups=None, n_features=20, feature_symbols=None,
                  key=None, fontsize=8, ncols=4, sharey=True,
                  show=None, save=None, ax=None, **kwds):
    sc.pl.rank_genes_groups(adata, groups=None, n_genes=n_features, gene_symbols=feature_symbols, key=None,
                            fontsize=8, ncols=4, sharey=True, show=None, save=None, ax=None, **kwds)


