import platform
if platform.system() != "Windows":
    #Note pysam doesn't support Windows
    import numpy as np
    import anndata as ad
    import pandas as pd
    import pysam
    from scipy.sparse import lil_matrix
    from tqdm import tqdm

    import scipy
    from .count import count


    def bld_mtx_fly(tsv_file, annotation, csv_file=None, genome=None, save=False):
        """
        Building count matrix on the fly.
        Expected running time for 10k cells X 100k features on a personal computer ~65min
        Does not count pcr duplicate.
        A tbi file with the same filename as the provided tsv_file must be present in the same directory as the tsv file
        Note that this function is not available on the Windows operating system.

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
        mtx = lil_matrix((nb_barcodes, len(window_list)), dtype=np.float32)
        for i, tmp_feat in enumerate(tqdm(window_list)):
            for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
                mtx[dict_barcodes[str(row).split('\t')[-2]], i] += 1

        print('building AnnData object')
        mtx = ad.AnnData(mtx.tocsr(),
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


def peak_mtx(fragments_file, peak_file, valid_bcs, normalized_peak_size=None, fast=False):

    names = ["chr", "start", "stop"]
    # cellranger peak files contain comments. So added comment parameter to skip them.
    features = pd.read_csv(peak_file, sep="\t", header=None, usecols=[0, 1, 2], names=names, comment='#')
    features.index = features.apply(lambda row: "_".join([str(val) for val in row]), axis=1)

    if normalized_peak_size:
        extension = int(np.ceil(normalized_peak_size / 2))
        start = round(features["start"] + (features["stop"] - features["start"]) / 2).astype(int) - extension
        stop = round(features["start"] + (features["stop"] - features["start"]) / 2).astype(int) + extension
        features["start"] = start
        features["stop"] = stop
        features["start"].clip(lower=0, inplace=True)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    ct_mtx = count(fragments_file, features.values.tolist(), valid_bcs, fast)

    if fast:
        X = scipy.sparse.csr_matrix(ct_mtx)
    else:
        X = ct_mtx.tocsr()

    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata


def gene_activity_mtx(fragments_file, gtf_file, valid_bcs, upstream=2000, downstream=0, source=None, gene_type=None, fast=False):

    names = ["chr", "source", "type", "start", "stop", "score", "strand", "frame", "attribute"]
    features = pd.read_csv(gtf_file, sep="\t", header=None, comment="#", names=names)

    features = features[features.type == "gene"]

    if source:
        features = features[features.source == source]

    features["gene_id"] = [attr.replace("gene_id", "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_id")]
    features["gene_name"] = [attr.replace("gene_name", "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_name")]
    features["gene_type"] = [attr.replace("gene_type", "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_type")]

    if gene_type:
        features = features[[feature in gene_type for feature in features.gene_type]]

    features.index = features.gene_id

    features["start"] = features.start - upstream
    features["start"].clip(lower=0, inplace=True)
    features["stop"] = features.stop + downstream

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    features = features[["gene_name", "gene_id", "gene_type", "chr", "start", "stop", "strand", "source"]]

    ct_mtx = count(fragments_file, features[["chr", "start", "stop"]].values.tolist(), valid_bcs, fast)

    if fast:
        X = scipy.sparse.csr_matrix(ct_mtx)
    else:
        X = ct_mtx.tocsr()

    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata


def window_mtx(fragments_file, valid_bcs, window_size=5000, species="human", fast=False):

    features = epi.ct.make_windows(window_size, chromosomes=species)

    features = [["chr{}".format(chrom), *window[:-1]] for chrom, windows in features.items() for window in windows]

    features = pd.DataFrame(features, columns=["chr", "start", "stop"])
    features.index = features.apply(lambda row: "_".join([str(val) for val in row]), axis=1)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    ct_mtx = count(fragments_file, features.values.tolist(), valid_bcs, fast)

    if fast:
        X = scipy.sparse.csr_matrix(ct_mtx)
    else:
        X = ct_mtx.tocsr()

    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata
