import anndata as ad
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import NMF
from sklearn.decomposition import LatentDirichletAllocation
from sklearn.datasets import make_multilabel_classification

import warnings
from warnings import warn


def tfidf(adata, norm='l1', layer_tfidf_key='tf-idf'):
    """
    TF-IDF --> Term Frequency Inverse Document Frequency
    scikit-learn implementation
    
    Parameters
    ----------
    
    adata : AnnData object
    
    norm : {'l1', 'l2'}, default='l2'
    Each output row will have unit norm, either:
    * 'l2': Sum of squares of vector elements is 1. The cosine
    similarity between two vectors is their dot product when l2 norm has
    been applied.
    * 'l1': Sum of absolute values of vector elements is 1.
    
    layer_tfidf_key : name of the tf-idf layer of the AnnData.
    
    """
    
    tfidf = TfidfTransformer(norm=norm)
    tfidf.fit(adata.X)
    tf_idf_matrix = tfidf.transform(adata.X)
    adata.layers[layer_tfidf_key] = tf_idf_matrix
    
def lsi(adata,
        n_components=50,
        tf_idf_layer='tf-idf',
        algorithm='arpack',
        n_iter=10,
        random_state=42):
    """
    Latent Semantinc Index -- > TF-IDF + TruncatedSVD.
    scikit-learn implementation
    
    Parameters
    ----------
    
    adata :
    
    n_components :
    
    tf_idf_layer :
    
    algorithm : arpack
    
    n_iter :
    
    random_state :
    
    """
    ## part 1 - TF-IDF
    if tf_idf_layer=='tf-idf' and (tf_idf_layer not in adata.layers.keys()):
        #add warning that 'tf-idf' was not computer and is getting computed with norm='L1' now
        warnings.warn("'tf-idf' layer doesn't exist and is now computed with norm='l1'")
        tf_idf_layer=None
        
    elif tf_idf_layer==None:
        adata.layers['tf-idf'] = tfidf(adata, norm='l1', layer_tfidf_key='tf-idf')
        tf_idf_layer='tf-idf'
    else:
        svd_model = adata.layers[tf_idf_layer]
        
    ## part 2 - decomposition
    svd_model = TruncatedSVD(n_components=n_components, 
                         algorithm=algorithm,
                         n_iter=n_iter, random_state=random_state)

    adata.obsm['X_lsi'] = svd_model.fit_transform(adata.layers[tf_idf_layer])



def nmf(adata, n_components=50, init='random', random_state=0):
	"""
	"""
	#start = time()
	model = NMF(n_components=n_components, init=init, random_state=random_state)
	adata.obsm['X_nmf'] = model.fit_transform(adata.X)
	adata.varm['NMF_components'] = model.components_.transpose()
	#end = time()
	#print(end-start)


    
def fa(adata, n_components=50, random_state=0):
	"""
	"""
    #start = time()
    transformer = FactorAnalysis(n_components=n_components, random_state=random_state)
    print('dense array... ', time()-start)
    adata.obsm['X_fa'] = transformer.fit_transform(adata.X.toarray())
    adata.varm['fa_components'] = transformer.components_.transpose()
    #end = time()
    #print(end-start)


def lda():
	"""
	"""
	#start = time()
	lda = LatentDirichletAllocation(n_components=20, random_state=0)
	adata.obsm['X_lda'] = lda.fit_transform(adata.X)
	#end = time()
	#print(end-start)





