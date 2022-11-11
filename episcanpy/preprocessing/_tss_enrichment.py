import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad


import multiprocessing
import pandas as pd

import platform
if platform.system() != "Windows":
    import pysam

def get_tss(gtf, source=None, feature=None, protein_coding_only=False):

    gtf_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'attribute', 'other']

    gtf = pd.read_csv(gtf,
        sep='\t',
        names=gtf_column_names,
        comment="#",
        header=None,
        index_col=False)

    if source is not None:
        gtf = gtf[gtf.source == source]

    if feature is not None:
        gtf = gtf[gtf.feature == feature]

    if protein_coding_only:
        gtf = gtf[gtf.other.str.contains("protein_coding")]

    tss_info = {"chromosome": [], "TSS_pos": [], "strand": []}

    for row in gtf.values:
        tss_info["chromosome"].append(row[0])
        tss_info["strand"].append(row[6])
        tss_info["TSS_pos"].append(row[3] if row[6] == "+" else row[4])

    return pd.DataFrame(data=tss_info)

# function that is executed by the processes (calculates coverage profiles for a subset of TSS)
def calc_coverage(preprocessed_input_chunk, mtx, bc_mapping, fragments_file, distance_to_tss, method, queue):

    # open fragments file
    tbx = pysam.TabixFile(fragments_file)

    # iterate over TSS subset
    for chromosome, tss_pos, strand in preprocessed_input_chunk:

        # calculate region boundaries
        region_start = tss_pos - distance_to_tss
        region_end = tss_pos + distance_to_tss + 1

        # iterate over reads overlapping the region
        for read in tbx.fetch(reference=chromosome, start=region_start, end=region_end, parser=pysam.asTuple()):

            barcode = read[3]

            # check if barcode is valid
            if barcode in bc_mapping:

                read_start = int(read[1])
                read_end = int(read[2])

                index_start = read_start - tss_pos + distance_to_tss
                index_end = read_end - tss_pos + distance_to_tss

                # check if index is in bounds
                if index_start < 0:
                    if method == "standard":
                        continue
                    elif method == "complete_coverage":
                        index_start = 0

                # check if index is in bounds
                if index_end >= mtx.shape[1]:
                    if method == "standard":
                        continue
                    elif method == "complete_coverage":
                        index_end = mtx.shape[1]

                # handle strand location
                if strand == "-":
                    tmp = mtx[bc_mapping[barcode]]
                    tmp = np.flip(tmp)

                    if method == "standard":
                        tmp[[index_start, index_end]] = tmp[[index_start, index_end]] + 1

                    elif method == "complete_coverage":
                        tmp[index_start:index_end] = tmp[index_start:index_end] + 1

                    mtx[bc_mapping[barcode]] = np.flip(tmp)

                else:
                    if method == "standard":
                        mtx[bc_mapping[barcode]][[index_start, index_end]] = mtx[bc_mapping[barcode]][[index_start, index_end]] + 1

                    elif method == "complete_coverage":
                        mtx[bc_mapping[barcode]][index_start:index_end] = mtx[bc_mapping[barcode]][index_start:index_end] + 1

    # put partial result in the queue
    queue.put(mtx)

def tss_enrichment(adata,
                   gtf,
                   fragments,
                   method="standard",
                   score="avg_score_of_center_region",
                   distance_to_tss=1000,
                   bp_per_flank=100,
                   n_jobs=1):
    """
    Computes the TSS profile for every observation
    The function offer two methods to compute the TSS enrichment coverage: 'standard' or 'complete_coverage'
    'standard' only consider the start and end point of the read.
    'complete_coverage' consider the entire read to calculate the coverage.

    Args:
    -----
        adata: AnnData
        gtf: path to GTF file
        fragments: path to fragments file
        method: method to use for computing TSS enrichment.
        score: value that is used as TSS enrichment score for individual observations
        distance_to_tss: distance to TSS
        bp_per_flank: number of base pair per flank to use for normalization
        n_jobs: number of jobs to use for computation. -1 means using all processors

    Returns:
    --------
        None
    """

    if method not in {"standard", "complete_coverage"}:
        raise ValueError("value for method argument must be 'standard' or 'complete_coverage'")

    if score not in {"avg_score_of_center_region", "score_at_zero"}:
        raise ValueError("value for score argument must be 'avg_score_of_center_region' or 'score_at_zero'")


    tss = get_tss(gtf, source="HAVANA", feature="gene", protein_coding_only=True)

    # if n_jobs == -1 use all available cores
    if n_jobs == -1:
        n_processes = multiprocessing.cpu_count()
    else:
        n_processes = n_jobs

    # create a dict that maps the barcode to a row number of the matrix
    bc_mapping = {bc: i for i, bc in enumerate(adata.obs_names)}

    # initialize matrix
    coverage_mtx = np.zeros((len(bc_mapping), distance_to_tss * 2 + 1))

    # create input chunks for the processes
    slice_size = int(np.ceil(tss.shape[0] / n_processes))
    max_idx = tss.shape[0]
    chunks = []

    for i in range(n_processes):

        # create slice indices
        if i == n_processes - 1:
            idxs = (slice_size * i, max_idx)
        else:
            idxs = (slice_size * i, slice_size * (i + 1))

        # preprocess TSS chunks
        chunk = []
        for chromosome, tss_pos, strand in zip(tss.chromosome[idxs[0]:idxs[1]],
                                               tss.TSS_pos[idxs[0]:idxs[1]],
                                               tss.strand[idxs[0]:idxs[1]]):

            region_start = tss_pos - distance_to_tss

            # skip mitochondrial TSS
            if chromosome == "chrM":
                continue
            # skip TSS too close to the start of the chromosome
            elif region_start < 0:
                continue
            # add valid TSS to the current chunk
            else:
                chunk.append((chromosome, tss_pos, strand))

        # add current chunk to the list of chunks
        chunks.append(chunk)


    # create queue
    q = multiprocessing.Queue()

    # spawn processes
    processes = []
    for i in range(n_processes):
        p = multiprocessing.Process(target=calc_coverage,
                                    args=(chunks[i], coverage_mtx.copy(), bc_mapping, fragments, distance_to_tss, method, q))
        p.start()
        processes.append(p)

    # remove and return items from queue
    results = []
    for _ in range(len(processes)):
        tmp = q.get()
        results.append(tmp)

    # block until all processes are terminated
    for p in processes:
        p.join()

    # put the partial results together
    for mtx in results:
        coverage_mtx += mtx

    # calculate the index of the center
    center = np.floor(coverage_mtx.shape[1] / 2)  # 0-index
    # create x-values
    x = [i - center for i in range(coverage_mtx.shape[1])]

    # create flank matrix
    coverage_mtx_flanks = np.concatenate((coverage_mtx[:, :bp_per_flank], coverage_mtx[:, -bp_per_flank:]), axis=1)

    # calculate mean coverage per barcode
    flank_means = coverage_mtx_flanks.mean(axis=1)
    flank_means = flank_means[:, None]

    # replace NAN and 0 flank means with population flank mean to avoid potential division by zero
    flank_means[np.isnan(flank_means)] = 0
    flank_means[flank_means == 0] = flank_means.mean()

    # calculate fold change of coverage over the mean coverage (flanks) per position per barcode
    fold_change_mtx = coverage_mtx / flank_means

    # calculate mean fold change per position
    mean_fold_change = fold_change_mtx.mean(axis=0)

    # TSS score equal the mean fold change at TSS position + (distance_to_tss / 2) positions before and after
    if score == "avg_score_of_center_region":
        tss_scores = fold_change_mtx[:, int(distance_to_tss / 2):(3 * int(distance_to_tss / 2) + 1)].mean(axis=1)
    # TSS score equal the fold change at TSS position
    elif score == "score_at_zero":
        tss_scores = fold_change_mtx[:, int(center)]

    mean_tss_score = tss_scores.mean()

    adata.obs["tss_enrichment_score"] = tss_scores

    adata.uns["tss_enrichment"] = {}
    adata.uns["tss_enrichment"]["x"] = x
    adata.uns["tss_enrichment"]["fold_change_mtx"] = fold_change_mtx
    adata.uns["tss_enrichment"]["mean_fold_change"] = mean_fold_change
    adata.uns["tss_enrichment"]["mean_tss_score"] = mean_tss_score


def tss_enrichment(adata,
                    group_by=None,
                    show_n=False,
                    max_cols=5,
                    figsize=None,
                    save=None):
    """
    Plots the mean TSS enrichment profile

    Args:
        adata: AnnData
        group_by: key of adata.obs that is used as grouping variable
        show_n: whether or not to show the number of observations
        max_cols: maximum number of columns
        figsize: size of the figure
        save: if True or str, save the figure. str is appended to the default filename. infer filetype if ending on {'.pdf', '.png', '.svg'}

    Returns:
        None
    """

    if group_by is None:

        if figsize is None:
            figsize = (6, 5)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        ax.plot(adata.uns["tss_enrichment"]["x"], adata.uns["tss_enrichment"]["mean_fold_change"], linewidth=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.grid(False)
        ax.set_xlabel("Distance from TSS (bp)")
        ax.set_ylabel("Mean TSS enrichment score")

    else:

        groups = np.sort(adata.obs[group_by].unique())
        n_groups = groups.size

        nrows = int(np.ceil(n_groups / max_cols))
        ncols = int(np.ceil(n_groups if n_groups <= max_cols else max_cols))

        if figsize is None:
            figsize = (ncols * 5 * 1.2, nrows * 5)

        fig, axs = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols, sharey=True, squeeze=False)
        axs = axs.flatten()

        for ax in axs[n_groups:]:
            ax.remove()

        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray",
                  "tab:olive", "tab:cyan"]

        for i, group in enumerate(groups):

            fold_change_mtx_part = adata.uns["tss_enrichment"]["fold_change_mtx"][adata.obs[group_by].str.fullmatch(group)]
            fold_change_part_mean = fold_change_mtx_part.mean(axis=0)

            ax = axs[i]

            ax.plot(adata.uns["tss_enrichment"]["x"], fold_change_part_mean, linewidth=0.5, c=colors[i] if n_groups <= 10 else colors[0])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.grid(False)
            ax.set_xlabel("Distance from TSS (bp)")
            ax.set_ylabel("Mean TSS enrichment score")
            ax.set_title(group, loc="center")

            if show_n:
                ax.annotate(text="n = {}".format(fold_change_mtx_part.shape[0]),
                            xy=(adata.uns["tss_enrichment"]["x"][0] + 0.05 * fold_change_part_mean.shape[0], 0.9 * fold_change_part_mean.max()),
                            xytext=(0, 0),
                            textcoords="offset points",
                            fontsize=12,
                            ha="left",
                            va="center")

    fig.suptitle("TSS enrichment")
    plt.tight_layout()

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "tss_enrichment.png"

        plt.savefig(filename, dpi=300)


def tss_enrichment_score(adata, figsize=(4, 6), save=None):
    """
    Plots a violin plot of the individual TSS enrichment scores

    Args:
        adata: AnnData
        figsize: size of the figure
        save: if True or str, save the figure. str is appended to the default filename. infer filetype if ending on {'.pdf', '.png', '.svg'}

    Returns:
        None
    """

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    sns.violinplot(y=adata.obs.tss_enrichment_score, inner=None, cut=0, ax=ax)
    sns.stripplot(y=adata.obs.tss_enrichment_score, size=1, jitter=0.3, color="black", ax=ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(False)
    ax.set_ylabel("TSS enrichment score")
    ax.set_title("TSS enrichment", loc="center")

    sns.despine()

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "tss_enrichment_score.png"

        plt.savefig(filename, dpi=300)



# splitting the cells by high and low enrichment score
def filter_enrichment_score(adata, score_threshold):
    """
    Cluster the cells in 2 categories (low and high tss enrichment) based on the mean tss enrichment score
    obtained with epi.pp.tss_enrichment()

    score <= score_threshold: --> low enrichment
    score > score thrshold --> high enrichment
    """
    annot = []
    for score in adata.obs['tss_enrichment_score']:
        if score <= score_threshold:
            annot.append('low enrichment')
        else:
            annot.append('high enrichment')
    adata.obs['tss_enrichment_split'] = annot
