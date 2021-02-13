from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import NMF
from sklearn.decomposition import TruncatedSVD
#from scanpy.plotting._tools.scatterplots import embedding
#from time import time


def nmf(adata, n_components=50, init='random', random_state=0):
    #start = time()
    model = NMF(n_components=n_components, init=init, random_state=random_state)
    adata.obsm['X_nmf'] = model.fit_transform(adata.X)
    adata.varm['NMF_components'] = model.components_.transpose()
    #end = time()
    #print(end-start)
    
def fa(adata, n_components=50, random_state=0):
    #start = time()
    transformer = FactorAnalysis(n_components=n_components, random_state=random_state)
    print('dense array... ', time()-start)
    adata.obsm['X_fa'] = transformer.fit_transform(adata.X.toarray())
    adata.varm['fa_components'] = transformer.components_.transpose()
    #end = time()
    #print(end-start)

def lsi(adata, n_components=50, n_iter=7, random_state=0):
    #start = time()
    svd = TruncatedSVD(n_components=n_components, n_iter=n_iter, random_state=random_state)
    output_X = svd.fit_transform(adata.X)
    adata.obsm['X_lsi'] = output_X
    adata.uns['lsi'] = {}
    adata.uns['lsi']['params'] = svd.get_params()
    adata.uns['lsi']['singular_values'] = svd.singular_values_
    adata.uns['lsi']['explained_variance_ratio'] = svd.explained_variance_ratio_
    #end = time()
    #print(end-start)