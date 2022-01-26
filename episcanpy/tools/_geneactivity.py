import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix
from intervaltree import Interval, IntervalTree
from collections import defaultdict


def _peaks_to_IntervalTree(adata):
    """
    Return a dict chr => Interval Tree
    Because the IntervalTree Lib is start inclusive end exclusive, I added +1 to end value.
    :return:dict[str] => IntervalTree
    """

    d_peaks = defaultdict(dict)
    for i, v in enumerate(adata.var_names.tolist()):

        spt = v.split("_")
        chrom, start, end = spt[0], int(spt[1]), int(spt[2])

        if chrom not in d_peaks:
            d_peaks[chrom] = IntervalTree()
        d_peaks[chrom].addi(start, end + 1, i)

    return d_peaks


def _test_genes_against_peaks(genes, peaks_intertree, matrix):

    ps = peaks_intertree[genes[0]: genes[1]]
    if ps :
        ps = list(ps)
        indice = [p[2] for p in ps]  # way faster to pull everything at once
        return (np.sum(matrix[:, indice].todense(), axis=1), genes[2])


def _read_gtf(gtf_file,
             upstream,
             feature_type,
             annotation):
    """
    Read a .gtf file return  dico[chromosome] => IntervalTree
    gtf_file: str || pathlib.Path
    upstream: int
    feature_type: str
    annotation: str

    """

    gtf = defaultdict(list)
    with open(gtf_file) as f:
        for line in f:
            if not line.startswith("#"):  # comment
                line = line.rstrip('\n').split('\t')

                if len(line) >= 8 and line[1] == annotation and line[2] == feature_type:

                    # negative strand
                    if line[6] == '-':
                        start, end = int(line[3]), int(line[4]) + upstream + 1
                    else:
                        start, end = int(line[3]) - upstream, int(line[4]) + 1

                    features = line[-1].rstrip(";").split(';')
                    feat_dict = {}
                    for y in features:
                        y = y.strip().split()
                        feat_dict[y[0].replace("'", "").replace('"', "")] = y[1].replace("'", "").replace('"', "")

                    gtf[line[0]].append([start, end, feat_dict])

    return gtf


def geneactivity(adata,
                 gtf_file,
                 key_added='gene',
                 upstream=5000,
                 feature_type='gene',
                 annotation='HAVANA',
                 layer_name='geneactivity',
                 raw=False,
                 copy=True):
    """

    Build an AnnData object containing the number of open features
    (windows, peaks, etc) overlapping genes (gene bodies + 5kb upstream of the TSS).
    It is possible to extend the distance from the TSS with the upstream parameter.

    Rather than using multiple time the same gene, it is possible to specify which genome annotation is desired using the parameter annotation.
    as the GTF files can often contain multiple annotations for genes (HAVANA, ENSEMBL, etc.).

    Alternatively, if you want to obtain the gene activity at something else than genes, like transcipts. It is possible as well.
    The feature_type can be specified.

    TSS = Transcription Starting Site

    INPUT
    -----

    adata : input AnnData
    gtf_file : input gtf file name + path
    key_added : unused / to save the geneactivity matrix as an adata.uns object if
    upstream :
    featyre_type : transcripts or genes
    annotation :
    layer_name : unused
    raw :
    copy : unused


    OUTPUT
    ------


    """
    gtf = _read_gtf(gtf_file=gtf_file,
                    upstream=upstream,
                    feature_type=feature_type,
                    annotation=annotation)

    dico_peaks = _peaks_to_IntervalTree(adata=adata)
    if raw:
        adata_raw = adata.raw.X.copy()
    else:
        adata_raw = adata.X.copy()

    m = csc_matrix(adata_raw)
    # could be easily adapt to multiple core. eg. st of cells of by chromosome or set of peaks
    genes_index = []
    genes = []
    for chrom in gtf:

        if chrom in dico_peaks:

            peaks_IT = dico_peaks[chrom]

            for gene_current in gtf[chrom]:

                r = _test_genes_against_peaks(gene_current,
                                                   peaks_IT,
                                                   m)
                if r :
                    genes_index.append(r[1])
                    genes.append(r[0])

    gene_activity = np.concatenate(genes, axis=-1)

    meta_data = set()

    for gene_annotation in genes_index:

        for annotation in gene_annotation.keys():
            meta_data.add(annotation)

    metadata_genes = {}
    for annot in meta_data:
        metadata_genes[annot] = []

    for gene_annotation in genes_index:
        for key in metadata_genes.keys():
            metadata_genes[key].append(gene_annotation.get(key, "NA"))


    for line in genes_index:
            dico_line = {}
            for element in line:
                if ' "' in element:
                    dico_line[element.rstrip('"').lstrip(" ").split(' "')[0]] = element.rstrip('"').lstrip(" ").split(' "')[1]

            for key in metadata_genes.keys():
                if key in dico_line.keys():
                    metadata_genes[key].append(dico_line[key])
                else:
                    metadata_genes[key].append('NA')

    # manage index
    dataframe_genes = pd.DataFrame.from_dict(metadata_genes)
    index_col = []
    if feature_type == 'transcript':
        if "transcript_name" in dataframe_genes.columns:
            index_col = ["transcript_name"]
        elif "transcript_id" in dataframe_genes.columns:
            index_col = ["transcript_id"]
    else:
        if "gene_name" in dataframe_genes.columns:
            index_col = ["gene_name"]
        elif "gene_id" in dataframe_genes.columns:
            index_col = ["gene_id"]

    if index_col:
        dataframe_genes.index = pd.Index(dataframe_genes[index_col[0]])
        dataframe_genes = dataframe_genes.drop(index_col, axis=1)

    gene_adata = ad.AnnData(gene_activity, var=dataframe_genes, obs=adata_raw.obs)
    gene_adata.uns = adata.uns.copy()
    gene_adata.obsm = adata.obsm.copy()
    gene_adata.obsp = adata.obsp.copy()
    return(gene_adata)
