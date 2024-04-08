import platform
if platform.system() != "Windows":
    #Note pysam doesn't support Windows
    import numpy as np
    import anndata as ad
    import pandas as pd
    import pysam
    from scipy.sparse import lil_matrix, csr_matrix, dok_matrix
    from intervaltree import Interval, IntervalTree
    from tqdm import tqdm
    import gzip
    import multiprocessing
    from ._features import make_windows


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



def get_barcodes(fragments,
                 comment="#"):

    if fragments.endswith(".gz"):
        fh = gzip.open(fragments, mode="rt")
    else:
        fh = open(fragments, mode="r")

    try:

        barcodes = set()

        check_for_comments = True
        use_strip = None

        # iterate over lines and count
        for i, line in enumerate(fh):

            # only check for comments at start - performance reasons
            if check_for_comments:
                if line.startswith(comment):
                    continue
                else:
                    check_for_comments = False

            # only strip if necessary - performance reasons
            if use_strip is None:
                line_split = line.strip().split("\t")
                if len(line_split) == 4:
                    use_strip = True
                else:
                    use_strip = False
            elif not use_strip:
                line_split = line.split("\t")
            else:
                line_split = line.strip().split("\t")

            # extract barcode
            bc = line_split[3]
            barcodes.add(bc)

    finally:
        fh.close()

    n_lines = i + 1

    return list(barcodes), n_lines



def count_lines(fragments):

    if fragments.endswith(".gz"):
        fh = gzip.open(fragments, mode="rt")
    else:
        fh = open(fragments, mode="r")

    try:
    
        # iterate over lines and count
        for i, line in enumerate(fh):
            pass

    finally:
        fh.close()

    n_lines = i + 1

    return n_lines



def count_fragments(fragments,
                    genomic_tree,
                    bc_mapping,
                    mode,
                    shape,
                    n_lines,
                    n_jobs):

    # initialize count matrix
    ct_mtx = csr_matrix(shape)

    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    # determine intervals for parallel processing
    lines_per_job = n_lines // n_jobs
    intervals = [(i*lines_per_job, (i+1)*lines_per_job) for i in range(n_jobs)]
    intervals[-1] = (intervals[-1][0], n_lines)

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = [CountProcess(fragments, genomic_tree, bc_mapping, mode, shape, intervals[i], q, sem) for i in range(n_jobs)]

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    # sum up partial count matrices
    while not q.empty():
        ct_mtx += q.get()

    return ct_mtx



class CountProcess(multiprocessing.Process):

    def __init__(self, fragments, genomic_tree, bc_mapping, mode, shape, interval, q, sem):

        super().__init__()

        self.fragments = fragments
        self.genomic_tree = genomic_tree
        self.bc_mapping = bc_mapping

        self.mode = mode
        self.shape = shape

        self.interval = interval

        self.q = q
        self.sem = sem


    def run(self):

        # initialize count matrix (DOK format for efficiency in incremental construction)
        ct_mtx = dok_matrix(self.shape)

        if self.fragments.endswith(".gz"):
            fh = gzip.open(self.fragments, mode="rt")
        else:
            fh = open(self.fragments, mode="r")

        try:

            for i, f in enumerate(fh):

                # skip lines outside of interval
                if i < self.interval[0]:
                    continue
                elif i >= self.interval[1]:
                    break
                    
                # skip comments
                if f.startswith("#"):
                    continue

                f = f.split("\t")

                bc = f[3]

                # skip invalid barcodes
                if bc not in self.bc_mapping:
                    continue

                chrom = f[0]
                start = int(f[1])
                stop = int(f[2])

                try:

                    # count based on transposition events
                    if self.mode == "transposition":

                        for pos in [start, stop]:
                            res = self.genomic_tree[chrom].at(pos)
                            feature_indices = [interval.data for interval in res]
                            for feature_idx in feature_indices:
                                ct_mtx[self.bc_mapping[bc], feature_idx] += 1

                    # count based on overlaps
                    elif self.mode == "overlap":
                        
                        res = self.genomic_tree[chrom].overlap(start, stop)
                        feature_indices = [interval.data for interval in res]
                        for feature_idx in feature_indices:
                            ct_mtx[self.bc_mapping[bc], feature_idx] += 1

                # skip unannotated chromosomes
                except KeyError:
                    continue

        finally:
            fh.close()

        # change to CSR format for efficiency in arithmetic operations and storage
        self.q.put(csr_matrix(ct_mtx))
        self.sem.release()



def peak_mtx(fragments_file,
             peak_file,
             valid_bcs=None,
             normalized_peak_size=None,
             mode="transposition",
             n_jobs=-1):
    """
    Generates a count matrix based on peaks.

    Args:
        fragments_file: path to fragments file
        peak_file: path to BED file
        valid_bcs: list of valid barcodes (optional)
        normalized_peak_size: if True peaks size will be normalized accordingly; default: None (no normalization)
        mode: determines what the counting is based on; either "transposition" or "overlap"; default: "transposition"
        n_jobs: number of parallel processes to start; default: -1 (all available cores)

    Returns:
        AnnData object
    """

    # load peaks
    names = ["chr", "start", "stop"]
    features = pd.read_csv(peak_file, sep="\t", header=None, usecols=[0, 1, 2], names=names, comment='#', dtype={"chr": str})

    # check if first row should be skipped (header)
    try:
        int(features.iloc[0].start)
    except ValueError:
        features = features[1:].copy()
        features["start"] = features.start.astype(int)
        features["stop"] = features.stop.astype(int)

    features.index = features.apply(lambda row: "_".join([str(val) for val in row]), axis=1)

    # normalize peak size
    if normalized_peak_size:
        extension = int(np.ceil(normalized_peak_size / 2))
        centers = ((features["start"] + features["stop"]) // 2).astype(int)
        start = centers - extension
        stop = centers + extension
        features["start"] = start
        features["stop"] = stop
        features["start"].clip(lower=0, inplace=True)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    # count lines to determine intevals for parallel processing; if valid_bcs=None, also get barcodes
    if valid_bcs is None:
        valid_bcs, n_lines = get_barcodes(fragments_file)
    else:
        n_lines = count_lines(fragments_file)

    # create a barcode mapping to associate barcodes with indices
    bc_mapping = {bc: i for i, bc in enumerate(valid_bcs)}

    mtx_shape = (len(valid_bcs), features.shape[0])  

    # create a genomic tree
    genomic_tree = {chrom: IntervalTree() for chrom in features.chr.unique()}
    for feature_idx, feature in enumerate(features.itertuples()):
        genomic_tree[feature.chr].add(Interval(feature.start, feature.stop, feature_idx))

    # create count matrix
    X = count_fragments(fragments_file, genomic_tree, bc_mapping, mode, mtx_shape, n_lines, n_jobs)

    # create AnnData object
    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata



def gene_activity_mtx(fragments_file,
                      gtf_file,
                      valid_bcs=None,
                      upstream=2000,
                      downstream=0,
                      source=None,
                      gene_type=None,
                      mode="transposition",
                      n_jobs=-1):
    """
    Generates a count matrix based on the openness of the gene bodies and promoter regions (gene activity).

    Args:
        fragments_file: path to fragments file
        gtf_file: path to GTF file
        valid_bcs: list of valid barcodes (optional)
        upstream: number of bp to consider upstream of TSS; default: 2000 bp
        downstream: number of bp to consider downstream of gene body; default: 0 bp
        source: filter for source of the feature; default: None (no filtering)
        gene_type: filter for gene type of the feature; default: None (no filtering)
        mode: determines what the counting is based on; either "transposition" or "overlap"; default: "transposition"
        n_jobs: number of parallel processes to start; default: -1 (all available cores)

    Returns:
        AnnData object
    """

    # load features from GTF file
    names = ["chr", "source", "type", "start", "stop", "score", "strand", "frame", "attribute"]
    features = pd.read_csv(gtf_file, sep="\t", header=None, comment="#", names=names, dtype={"chr": str})

    # filter for genes
    features = features[features.type == "gene"]

    # filter for source
    if source:
        features = features[features.source == source]

    # extract gene ID
    features["gene_id"] = [attr.replace("gene_id", "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_id")]

    # extract gene name
    possible_tags = ["gene_name", "gene"]
    for tag in possible_tags:
        tmp = [attr.replace(tag, "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith(f"{tag} ")]
        if tmp:
            break
    features["gene_name"] = tmp

    # extract gene type
    possible_tags = ["gene_type", "gene_biotype"]
    for tag in possible_tags:
        tmp = [attr.replace(tag, "").strip().strip("\"") for feature_attr in features.attribute for attr in feature_attr.split(";") if attr.strip().startswith(f"{tag} ")]
        if tmp:
            break
    features["gene_type"] = tmp

    # filter for gene type
    if gene_type:
        features = features[[feature in gene_type for feature in features.gene_type]]

    features.index = features.gene_id

    # adjust coordinates for upstream and downstream parameters
    features["start"] = features.apply(lambda row: row.start - upstream if row.strand == "+" else row.start - downstream, axis=1)
    features["start"].clip(lower=0, inplace=True)
    features["stop"] = features.apply(lambda row: row.stop + downstream if row.strand == "+" else row.stop + upstream, axis=1)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    features = features[["gene_name", "gene_id", "gene_type", "chr", "start", "stop", "strand", "source"]]

    # count lines to determine intevals for parallel processing; if valid_bcs=None, also get barcodes
    if valid_bcs is None:
        valid_bcs, n_lines = get_barcodes(fragments_file)
    else:
        n_lines = count_lines(fragments_file)

    # create a barcode mapping to associate barcodes with indices
    bc_mapping = {bc: i for i, bc in enumerate(valid_bcs)}
    
    mtx_shape = (len(valid_bcs), features.shape[0]) 

    # create a genomic tree
    genomic_tree = {chrom: IntervalTree() for chrom in features.chr.unique()}
    for feature_idx, feature in enumerate(features.itertuples()):
        genomic_tree[feature.chr].add(Interval(feature.start, feature.stop, feature_idx))

    # create count matrix
    X = count_fragments(fragments_file, genomic_tree, bc_mapping, mode, mtx_shape, n_lines, n_jobs)

    # create AnnData object
    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata



def window_mtx(fragments_file,
               valid_bcs=None,
               window_size=5000,
               species="human",
               mode="transposition",
               n_jobs=-1):
    """
    Generates a count matrix based on the openness of equally sized bins of the genome (windows).

    Args:
        fragments_file: path to fragments file
        valid_bcs: list of valid barcodes (optional)
        window_size: size of windows in bp; default: 5000 bp
        species: species to create the windows for (human or mouse); default: "human"; will be extended in the future
        mode: determines what the counting is based on; either "transposition" or "overlap"; default: "transposition"
        n_jobs: number of parallel processes to start; default: -1 (all available cores)

    Returns:
        AnnData object
    """

    # create windows
    features = make_windows(window_size, chromosomes=species)

    features = [["chr{}".format(chrom), *window[:-1]] for chrom, windows in features.items() for window in windows]

    # create feature DataFrame
    features = pd.DataFrame(features, columns=["chr", "start", "stop"])
    features["chr"] = features.chr.astype(str)

    features.index = features.apply(lambda row: "_".join([str(val) for val in row]), axis=1)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    # count lines to determine intevals for parallel processing; if valid_bcs=None, also get barcodes
    if valid_bcs is None:
        valid_bcs, n_lines = get_barcodes(fragments_file)
    else:
        n_lines = count_lines(fragments_file)

    # create a barcode mapping to associate barcodes with indices
    bc_mapping = {bc: i for i, bc in enumerate(valid_bcs)}
    
    mtx_shape = (len(valid_bcs), features.shape[0]) 

    # create a genomic tree
    genomic_tree = {chrom: IntervalTree() for chrom in features.chr.unique()}
    for feature_idx, feature in enumerate(features.itertuples()):
        genomic_tree[feature.chr].add(Interval(feature.start, feature.stop, feature_idx))

    # create count matrix
    X = count_fragments(fragments_file, genomic_tree, bc_mapping, mode, mtx_shape, n_lines, n_jobs)

    # create AnnData object
    adata = ad.AnnData(X, obs=pd.DataFrame(index=valid_bcs), var=features)

    return adata
