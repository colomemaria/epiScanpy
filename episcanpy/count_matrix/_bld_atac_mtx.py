import numpy as np
import anndata as ad
import pandas as pd
import pysam
from scipy.sparse import csc_matrix
from tqdm import tqdm


def bld_mtx_fly(tsv_file, annotation, csv_file=None, genome=None, save=False):
    """
    Building count matrix on the fly.
    Expected running time for 10k cells X 100k features on a personal computer ~65min
    Does not count pcr duplicate.
    A tbi file with the same filename as the provided tsv_file must be present in the same directory as the tsv file

    Parameters
    ----------

    tsv_file : name of the file containing the multiplexed reads (.tsv or .tsv.gz)

    annotation : loaded set of features to create the feature matrix from

    csv_file : default is None -

    genome : default is None - specify if you want to extract a specific genome assembly

    save : default is False - supply a file path as str to save generated AnnData object

    Output
    ------

    AnnData object (also saved as h5ad if save argument is specified)

    """

    print('loading barcodes')
    barcodes = sorted(pd.read_csv(tsv_file, sep='\t', header=None).loc[:, 3].unique().tolist())

    # barcodes
    nb_barcodes = len(barcodes)
    dict_barcodes = {barcodes[i]: i for i in range(0, len(barcodes))}

    # Load tabix
    tbx = pysam.TabixFile(tsv_file)

    # format annotations
    window_list = []

    if genome:
        for chrom in sorted(annotation.keys()):
            window_list += [["".join([genome, '_chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]
    else:
        for chrom in sorted(annotation.keys()):
            window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]

    print('building count matrix')
    mtx = []
    for tmp_feat in tqdm(window_list):
        vector = [0]*nb_barcodes
        for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
            vector[dict_barcodes[str(row).split('\t')[-2]]] += 1
        mtx.append(vector)

    print('building AnnData object')
    mtx = csc_matrix(np.array(mtx))
    mtx = ad.AnnData(mtx.transpose(),
                     obs=pd.DataFrame(index=barcodes),
                     var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in window_list]))

    if csv_file:
        print('filtering barcodes')
        df = pd.read_csv(csv_file)
        if genome == 'mm10':
            df_filtered = df[(df.is_mm10_cell_barcode == 1)]
        elif genome == 'hg19':
            df_filtered = df[(df.is_hg19_cell_barcode == 1)]
        else:
            df_filtered = df[(df.is_cell_barcode == 1)]

        barcodes = set(df_filtered.barcode.tolist())
        mtx = mtx[[i in barcodes for i in mtx.obs.index]].copy()

    if save:
        mtx.write(save)

    return mtx
