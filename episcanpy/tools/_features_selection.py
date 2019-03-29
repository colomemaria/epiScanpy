import scanpy.api as sc

def rank_features(adata, groupby, use_raw=True, groups='all', reference='rest', n_genes=100,
                     rankby_abs=False, key_added=None, copy=False, method='',
                     corr_method='benjamini-hochberg', **kwds):

    if (method!='') and (adata.uns['omic'] != 'methylation'):
        method='t-test_overestim_var'
    else:
        method='t-test'

    sc.tl.rank_genes_groups(adata, groupby, use_raw=True, groups='all', reference='rest', n_genes=100,
                     rankby_abs=False, key_added=None, copy=False, method='',
                     corr_method='benjamini-hochberg', **kwds)
