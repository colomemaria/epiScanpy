#from scanpy.api.tl import rank_genes_groups
import scanpy as sc
import warnings
from warnings import warn

def rank_features(adata, groupby, omic=None, use_raw=True, groups='all', reference='rest', n_features=100,
                     rankby_abs=False, key_added='rank_features_groups', copy=False, method='',
                     corr_method='benjamini-hochberg', **kwds):

    """
    It is a wrap-up function of scanpy sc.tl.rank_genes_groups function. For more information see 
    Scanpy documentations.

    Here, we optimised the default features for the different kind of omics epiScanpy is analysing.
    Parameters for rankby_abs, method and corr_method are fixed byt the omic (if not specified otherwise)

    In case you want to change all the settings, it is advise to directly use the Scanpy function 
    or to specify the omic parameter as False. 

    If the omic of the input AnnData object is not specified (or incorrect), you can add it in omic 
    (either 'RNA', 'ATAC' or 'methylation'). If the omic of the current matrix is not known by
    epiScanpy but you wnat to use settings of a known omic, specify which omic as a parameter. 

    """
    if omic == None:
    	if 'omic' in adata.uns.keys():
    		omic = adata.uns['omic']
    	else:
    		warn("""Attention: no omic specified. We used default settings of the original Scanpy function\n
    			When the parameters where not specified in input""")
    		omic = 'RNA'

    if (method!='') and (adata.uns['omic'] != 'methylation'):
        method='t-test_overestim_var'
    else:
        method='t-test'
    if omic == 'methylation':
    	if copy==False:
    		sc.tl.rank_genes_groups(adata=adata, groupby=groupby, use_raw=use_raw,
                groups=groups, reference=reference, n_genes=n_features,
    			rankby_abs=True, key_added=key_added, copy=False, method='t-test', corr_method='benjamini-hochberg')
    		return()
    	else:
    		adata2 = sc.tl.rank_genes_groups(adata=adata, groupby=groupby, use_raw=use_raw,
                groups=groups, reference=reference, n_genes=n_features,
    			rankby_abs=True, key_added=key_added, copy=True, method='t-test', corr_method='benjamini-hochberg')
    		return(adata2)
    else:
    	if copy==False:
    		sc.tl.rank_genes_groups(adata=adata, groupby=groupby, use_raw=use_raw,
                            groups=groups, reference=reference, n_genes=n_features,
                     		rankby_abs=rankby_abs, key_added=key_added, copy=False, method=method,
                     		corr_method=corr_method, **kwds)
    		return()
    	else:
    		adata2 = sc.tl.rank_genes_groups(adata=adata, groupby=groupby, use_raw=use_raw,
                            groups=groups, reference=reference, n_genes=n_features,
                     		rankby_abs=rankby_ab, key_added=key_added, copy=True, method=method,
                            corr_method=corr_method, **kwds)
    		return(adata2)


