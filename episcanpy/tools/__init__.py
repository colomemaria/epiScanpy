#from ..tools._impute_meth import readandimputematrix
from ..tools._features_selection import rank_features
from ..tools._silhouette import silhouette
from ..tools._recipe import lazy
from ..tools._cell_id import load_markers, identify_cluster
from ..tools._top_feature_genes import top_feature_genes, var_features_to_genes
#from ..tools._find_genes2 import find_genes_in_features
from ..tools._scanpy_fct import pca, diffmap, draw_graph, tsne, umap
#from ..tools._scanpy_fct import louvain, leiden
from ..tools._scanpy_fct import dendogram
from ..tools._find_genes import find_genes
from scanpy.tools import dpt

from ..tools._clustering import getNClusters
from ..tools._clustering import kmeans, hc
from ..tools._clustering import ARI, AMI, homogeneity
from scanpy.tools import louvain, leiden
from ..tools._geneactivity import geneactivity

from ..tools._impute_gene_methylation import imputation_feature
from ..tools._comparisons import transfer_obs, imputation

from ..preprocessing._decomposition import tfidf, lsi
from ..preprocessing._decomposition import nmf, fa