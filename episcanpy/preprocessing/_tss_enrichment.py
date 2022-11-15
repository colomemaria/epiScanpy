import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad

import gzip
import pandas as pd

import platform
if platform.system() != "Windows":
    import pysam


def get_tss(gtf,
            source=None,
            feature=None,
            protein_coding_only=False):

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

    tss_info = {"chr": [], "tss_pos": [], "strand": []}

    for row in gtf.values:
        tss_info["chr"].append(row[0])
        tss_info["strand"].append(row[6])
        tss_info["tss_pos"].append(row[3] if row[6] == "+" else row[4])

    return pd.DataFrame(data=tss_info)


def tss_enrichment(adata,
                   fragments,
                   gtf,
                   n=5000,
                   score="avg_score_of_center_region",
                   distance_to_tss=1000,
                   bp_per_flank=100,
                   comment="#"):
    """
    Computes the TSS profile for every observation.

    Args:
        adata: AnnData
        fragments: path to fragments file
        gtf: path to GTF file
        n: number of TSS to use for calculation
        score: value that is used as TSS enrichment score for individual observations
        distance_to_tss: distance to TSS
        bp_per_flank: number of base pair per flank to use for normalization
        comment: character that marks comments

    Returns:
        None
    """

    features = get_tss(gtf, source="HAVANA", feature="gene", protein_coding_only=True)

    features["start"] = features.tss_pos - distance_to_tss
    features["stop"] = features.tss_pos + distance_to_tss

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    features = features[["chr", "start", "stop", "tss_pos", "strand"]].values.tolist()

    if n:
        features = features[:n]

    check_for_comments = True
    use_strip = None

    feature_idx = 0

    chrom_mapping = dict()
    bc_mapping = {bc: i for i, bc in enumerate(adata.obs_names)}

    valid_bcs_set = set(bc_mapping.keys())

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    # initialize matrix
    count_mtx = np.zeros((len(bc_mapping), distance_to_tss * 2 + 1))

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    # mapping chromosome/contig names to numeric values
    c = 0
    for chrom in [feature[0] for feature in features]:
        if not chrom in chrom_mapping:
            chrom_mapping[chrom] = c
            c += 1

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    if fragments.endswith(".gz"):
        fh = gzip.open(fragments, mode="rt")
    else:
        fh = open(fragments, mode="r")

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    for line in fh:

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

        ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

        bc = line_split[3]
        if bc not in valid_bcs_set:
            continue

        ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

        chrom = line_split[0]

        # map chromosome/contig name to numeric value
        # skip if no features on chromosome/contig
        if chrom in chrom_mapping:
            chrom_int = chrom_mapping[chrom]
        else:
            continue

        # fragment on previous chromosome
        if chrom_int < chrom_mapping[features[feature_idx][0]]:
            continue
        # feature on previous chromosome
        elif chrom_int > chrom_mapping[features[feature_idx][0]]:
            while feature_idx < len(features) and chrom_int > chrom_mapping[features[feature_idx][0]]:
                feature_idx += 1
            if feature_idx >= len(features):
                break

        ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

        start = int(line_split[1])
        stop = int(line_split[2])

        # fragment in front of feature
        if stop < features[feature_idx][1]:
            continue

        # fragment behind feature
        elif start > features[feature_idx][2]:
            while feature_idx < len(features) and start > features[feature_idx][2] and chrom_int == chrom_mapping[features[feature_idx][0]]:
                feature_idx += 1

            # end of features
            if feature_idx >= len(features):
                break

            # fragment on previous chromosome
            elif chrom_int != chrom_mapping[features[feature_idx][0]]:
                continue

            # overlap
            elif not stop < features[feature_idx][1]:

                # full overlap
                if start >= features[feature_idx][1] and stop <= features[feature_idx][2]:

                    index_start = start - features[feature_idx][3] + distance_to_tss
                    index_stop = stop - features[feature_idx][3] + distance_to_tss

                    if features[feature_idx][4] == "-":
                        tmp = count_mtx[bc_mapping[bc]]
                        tmp = np.flip(tmp)
                        tmp[[index_start, index_stop]] = tmp[[index_start, index_stop]] + 1
                        count_mtx[bc_mapping[bc]] = np.flip(tmp)

                    else:
                        count_mtx[bc_mapping[bc]][[index_start, index_stop]] = count_mtx[bc_mapping[bc]][[index_start, index_stop]] + 1

                tmp_idx = feature_idx + 1
                is_overlapping = True
                on_same_chrom = True

                # follow-up
                while tmp_idx < len(features) and is_overlapping and on_same_chrom:
                    is_overlapping = not (stop < features[tmp_idx][1] or start > features[tmp_idx][2])
                    on_same_chrom = chrom_int == chrom_mapping[features[tmp_idx][0]]

                    if is_overlapping and on_same_chrom:
                        # full overlap
                        if start >= features[tmp_idx][1] and stop <= features[tmp_idx][2]:

                            index_start = start - features[tmp_idx][3] + distance_to_tss
                            index_stop = stop - features[tmp_idx][3] + distance_to_tss

                            if features[tmp_idx][4] == "-":
                                tmp = count_mtx[bc_mapping[bc]]
                                tmp = np.flip(tmp)
                                tmp[[index_start, index_stop]] = tmp[[index_start, index_stop]] + 1
                                count_mtx[bc_mapping[bc]] = np.flip(tmp)

                            else:
                                count_mtx[bc_mapping[bc]][[index_start, index_stop]] = count_mtx[bc_mapping[bc]][[index_start, index_stop]] + 1

                    tmp_idx += 1

        # overlap
        else:

            # full overlap
            if start >= features[feature_idx][1] and stop <= features[feature_idx][2]:

                index_start = start - features[feature_idx][3] + distance_to_tss
                index_stop = stop - features[feature_idx][3] + distance_to_tss

                if features[feature_idx][4] == "-":
                    tmp = count_mtx[bc_mapping[bc]]
                    tmp = np.flip(tmp)
                    tmp[[index_start, index_stop]] = tmp[[index_start, index_stop]] + 1
                    count_mtx[bc_mapping[bc]] = np.flip(tmp)

                else:
                    count_mtx[bc_mapping[bc]][[index_start, index_stop]] = count_mtx[bc_mapping[bc]][[index_start, index_stop]] + 1

            tmp_idx = feature_idx + 1
            is_overlapping = True
            on_same_chrom = True

            # follow-up
            while tmp_idx < len(features) and is_overlapping and on_same_chrom:
                is_overlapping = not (stop < features[tmp_idx][1] or start > features[tmp_idx][2])
                on_same_chrom = chrom_int == chrom_mapping[features[tmp_idx][0]]

                if is_overlapping and on_same_chrom:
                    # full overlap
                    if start >= features[tmp_idx][1] and stop <= features[tmp_idx][2]:

                        index_start = start - features[tmp_idx][3] + distance_to_tss
                        index_stop = stop - features[tmp_idx][3] + distance_to_tss

                        if features[tmp_idx][4] == "-":
                            tmp = count_mtx[bc_mapping[bc]]
                            tmp = np.flip(tmp)
                            tmp[[index_start, index_stop]] = tmp[[index_start, index_stop]] + 1
                            count_mtx[bc_mapping[bc]] = np.flip(tmp)

                        else:
                            count_mtx[bc_mapping[bc]][[index_start, index_stop]] = count_mtx[bc_mapping[bc]][[index_start, index_stop]] + 1

                tmp_idx += 1

    fh.close()

    # calculate the index of the center
    center = np.floor(count_mtx.shape[1] / 2)  # 0-index
    # create x-values
    x = [i - center for i in range(count_mtx.shape[1])]

    # create flank matrix
    count_mtx_flanks = np.concatenate((count_mtx[:, :bp_per_flank], count_mtx[:, -bp_per_flank:]), axis=1)

    # calculate mean coverage per barcode
    flank_means = count_mtx_flanks.mean(axis=1)
    flank_means = flank_means[:, None]

    # replace NAN and 0 flank means with population flank mean to avoid potential division by zero
    flank_means[np.isnan(flank_means)] = 0
    flank_means[flank_means == 0] = flank_means.mean()

    # calculate fold change of coverage over the mean coverage (flanks) per position per barcode
    fold_change_mtx = count_mtx / flank_means

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


def tss_enrichment_plot(adata,
                        group_by=None,
                        show_n=False,
                        max_cols=5,
                        figsize=None,
                        save=None):
    """
    Plots the mean TSS enrichment profile.

    Args:
        adata: AnnData
        group_by: key of .obs that is used as grouping variable
        show_n: whether show the number of observations
        max_cols: maximum number of columns
        figsize: size of the figure
        save: if True or str, save the figure. str represents entire path. filetype is inferred.

    Returns:
        None
    """

    if group_by is None:

        if figsize is None:
            figsize = (4, 4)

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
            figsize = (ncols * 4 * 1.2, nrows * 4)

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
