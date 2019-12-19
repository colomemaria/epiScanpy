from ..plotting._heatmap import overlap_heatmap, prct_overlap
#from ..plotting._features import plot_rank_features
from ..plotting._silhouette import silhouette, silhouette_tot
from ..plotting._histograms import cluster_composition

from ..plotting._scanpy_plotting import pca, tsne, umap, diffmap, draw_graph
from ..plotting._scanpy_plotting import pca_loadings, pca_overview, pca_variance_ratio
from ..plotting._scanpy_plotting import rank_feat_groups
from ..plotting._scanpy_plotting import rank_feat_groups_violin, rank_feat_groups_dotplot
from ..plotting._scanpy_plotting import rank_feat_groups_stacked_violin, rank_feat_groups_matrixplot
from ..plotting._scanpy_plotting import rank_feat_groups_heatmap, rank_feat_groups_tracksplot
from ..plotting._scanpy_plotting import dendrogram, correlation_matrix

import scanpy as sc
from scanpy.api.pl import scatter, violin, ranking, clustermap, stacked_violin, heatmap, dotplot, matrixplot, tracksplot

from ..preprocessing._quality_control import cal_var, variability_features

