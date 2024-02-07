import multiprocessing
import gzip

from collections import defaultdict

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import pandas as pd

from Bio import SeqIO

from ._motif_search import get_sequence




def get_features(adata, distance_to_center=100):

    features = defaultdict(list)

    for res in adata.uns["motif_search"]["results"]:

        region = adata.var.iloc[res["feature_idx"]]

        if res["pos"] > 0:
            motif_center = region.start + res["pos"] + np.floor(res["motif_length"] / 2)
        else:
            motif_center = region.stop + res["pos"] + np.ceil(res["motif_length"] / 2)

        features["chr"].append(region.chr)
        features["motif_center"].append(int(motif_center))
        features["strand"].append("+" if res["pos"] > 0 else "-")
        features["motif_idx"].append(res["motif_idx"])
        features["start"].append(int(motif_center) - distance_to_center)
        features["stop"].append(int(motif_center) + distance_to_center)

    features = pd.DataFrame(data=features)

    features.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)
    features.reset_index(inplace=True, drop=True)
    
    return features



def count_lines(filename):
    """Count the number of lines in a file."""

    if filename.endswith(".gz"):
        fh = gzip.open(filename, mode="rt")
    else:
        fh = open(filename, mode="r")

    for i, _ in enumerate(fh):
        pass

    fh.close()

    return i + 1



def cutsite_kmers(fragments,
                  genome,
                  adata,
                  n_lines,
                  kmer_size=6,
                  comment="#",
                  n_jobs=1):
    
    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    lines_per_job = n_lines // n_jobs
    intervals = [(i*lines_per_job, (i+1)*lines_per_job) for i in range(n_jobs)]
    intervals[-1] = (intervals[-1][0], n_lines)

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = [SeqReaderCutsites(fragments, genome, set(adata.obs_names), intervals[i], q, sem, kmer_size=kmer_size, comment=comment) for i in range(n_jobs)]

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    kmers = defaultdict(int)
    while not q.empty():
        for bc, count in q.get().items():
            kmers[bc] += count

    return kmers



class SeqReaderCutsites(multiprocessing.Process):

    def __init__(self, fragments, genome, valid_bcs, interval, q, sem, kmer_size=6, comment="#"):
        super().__init__()

        self.fragments = fragments
        self.genome = genome
        self.valid_bcs = valid_bcs

        self.kmer_size = kmer_size

        self.interval = interval
        self.comment = comment

        self.q = q
        self.sem = sem


    def run(self) -> None:
        
        chromosome_seq = None
        current_chromosome = None

        kmers = defaultdict(int)

        ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

        if self.fragments.endswith(".gz"):
            fh = gzip.open(self.fragments, mode="rt")
        else:
            fh = open(self.fragments, mode="r")

        ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

        for i, line in enumerate(fh):

            if i < self.interval[0]:
                continue
            elif i >= self.interval[1]:
                break

            if line.startswith(self.comment):
                continue

            line_split = line.strip().split("\t")

            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

            bc = line_split[3]
            if bc not in self.valid_bcs:
                continue

            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

            chrom = line_split[0]
            start = int(line_split[1])
            stop = int(line_split[2])

            # if chrom != current_chromosome:
            #     chromosome_seq = get_chrom_sequence(chrom, self.genome)
            #     current_chromosome = chrom

            # for pos in [start, stop]:
            #     kmer = chromosome_seq[ int(pos - np.floor(self.kmer_size/2)) : int(pos + np.ceil(self.kmer_size/2)) ]
            #     kmers[kmer] += 1

            for pos in [start, stop]:
                try:
                    kmer = get_sequence(chrom=chrom, start=int(pos - np.floor(self.kmer_size/2)), end=int(pos + np.ceil(self.kmer_size/2)), genome=self.genome)
                    kmers[kmer] += 1
                except KeyError:
                    continue

        self.q.put(kmers)

        self.sem.release()



def peak_bg_kmers(adata, kmer_size, genome):

    shift_left = int(np.floor(kmer_size / 2))
    shift_right = int(np.ceil(kmer_size / 2))

    kmers = defaultdict(int)
    for row in adata.var.to_dict("records"):

        try:
            seq = get_sequence(row["chr"], row["start"] - shift_left, row["stop"] + shift_right + 1, genome=genome)
        except KeyError:
            continue

        for i in range(shift_left, len(seq) - shift_right):
            window_start = i - shift_left
            window_stop = i + shift_right
            kmer = seq[window_start:window_stop]
            kmers[kmer] += 1

    return kmers



def genome_bg_kmers(genome, kmer_size, n_jobs=-1):

    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = []

    for chrom in genome.keys():

        p = KmerCounter(
            chrom,
            genome,
            kmer_size,
            sem=sem,
            q=q
        )

        producers.append(p)    

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    kmer_freqs = defaultdict(int)
    while not q.empty():
        kmer_freqs_partial = q.get()
        for kmer, count in kmer_freqs_partial.items():
            kmer_freqs[kmer] += count

    return kmer_freqs



class KmerCounter(multiprocessing.Process):

    def __init__(self, chrom, genome, kmer_size, sem, q):

        super().__init__()

        self.chrom = chrom
        self.genome = genome
        self.kmer_size = kmer_size

        self.sem = sem
        self.q = q


    def run(self):

        kmer_freqs_partial = defaultdict(int)

        seq = self.genome[self.chrom].seq

        for i in range(len(seq) - self.kmer_size):
            kmer = seq[i:i+self.kmer_size]
            kmer_freqs_partial[kmer] += 1

        self.q.put(kmer_freqs_partial)
        self.sem.release()



def build_tree_parallel(features, n_jobs=1):

    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = [ChromTree(features[features.chr == chrom], chrom, q, sem) for chrom in features.chr.unique()]

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    results = []
    while not q.empty():
        results.append(q.get())

    tree = dict()
    for subtree, chrom in results:
        tree[chrom] = dict(subtree)
        
    return tree



class ChromTree(multiprocessing.Process):

    def __init__(self, features, chrom, q, sem):
        super().__init__()

        self.features = features
        self.chrom = chrom

        self.q = q
        self.sem = sem


    def run(self) -> None:
        
        subtree = defaultdict(list)

        for feat_idx, feature in self.features.iterrows():

            for pos in range(feature.start, feature.stop + 1):
                subtree[pos].append(feat_idx)

        self.q.put((subtree, self.chrom))

        self.sem.release()



def tn5_activity(adata, fragments, features, groupby, pos_to_feat, genome, k, bias_correction_factors, n_lines, region_size=201, comment="#", n_jobs=1):

    group_mapping = {group: i for i, group in enumerate(adata.obs[groupby].cat.categories)}
    bc_mapping = {bc: group_mapping[group] for bc, group in zip(adata.obs_names, adata.obs[groupby])}

    n_motifs = len(features.motif_idx.unique())
    n_groups = len(adata.obs[groupby].unique())

    shape = (n_motifs, n_groups, region_size)
    ct_mtx_raw = np.zeros(shape)
    ct_mtx_corr = np.zeros(shape)

    features = features.values.tolist()

    
    n_jobs = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs

    lines_per_job = n_lines // n_jobs
    intervals = [(i*lines_per_job, (i+1)*lines_per_job) for i in range(n_jobs)]
    intervals[-1] = (intervals[-1][0], n_lines)

    sem = multiprocessing.Semaphore(n_jobs)

    manager = multiprocessing.Manager()
    q = manager.Queue()

    producers = [Tn5ActivityProcess(fragments, features, groupby, pos_to_feat, genome, k, bias_correction_factors, bc_mapping, shape, intervals[i], q, sem) for i in range(n_jobs)]

    for p in producers:
        sem.acquire()
        p.start()

    for p in producers:
        p.join()

    while not q.empty():
        ct_mtx_raw_partial, ct_mtx_corr_partial = q.get()
        ct_mtx_raw += ct_mtx_raw_partial
        ct_mtx_corr += ct_mtx_corr_partial
    

    adata.uns["footprinting"] = {}
    adata.uns["footprinting"]["group_mapping"] = {v: k for k, v in group_mapping.items()}
    adata.uns["footprinting"]["mtx_raw"] = ct_mtx_raw
    adata.uns["footprinting"]["mtx_corr"] = ct_mtx_corr



class Tn5ActivityProcess(multiprocessing.Process):

    def __init__(self, fragments, features, groupby, pos_to_feat, genome, k, bias_correction_factors, bc_mapping, shape, interval, q, sem):

        super().__init__()

        self.fragments = fragments
        self.features = features
        self.groupby = groupby
        self.pos_to_feat = pos_to_feat
        self.genome = genome
        self.k = k
        self.bias_correction_factors = bias_correction_factors
        self.bc_mapping = bc_mapping
        self.shape = shape
        self.interval = interval
        self.q = q
        self.sem = sem


    def run(self):

        ct_mtx_raw = np.zeros(self.shape)
        ct_mtx_corr = np.zeros(self.shape)
        
        if self.fragments.endswith(".gz"):
            fh = gzip.open(self.fragments, mode="rt")
        else:
            fh = open(self.fragments, mode="r")

        for i, f in enumerate(fh):

            if i < self.interval[0]:
                continue
            elif i >= self.interval[1]:
                break
                
            if f.startswith("#"):
                continue

            f = f.split("\t")

            bc = f[3]

            if bc not in self.bc_mapping:
                continue

            chrom = f[0]

            if chrom in self.pos_to_feat:

                start = int(f[1])
                stop = int(f[2])

                for pos in [start, stop]:

                    if pos in self.pos_to_feat[chrom]:

                        for feat_idx in self.pos_to_feat[chrom][pos]:

                            feature = self.features[feat_idx]

                            feature_start = feature[4]
                            feature_stop = feature[5]
                            feature_strand = feature[2]
                            feature_motif_idx = feature[3]

                            if feature_strand == "+":
                                rel_pos = pos - feature_start
                            else:
                                rel_pos = feature_stop - pos

                            kmer = get_sequence(chrom=chrom, start=int(pos - np.floor(self.k/2)), end=int(pos + np.ceil(self.k/2)), genome=self.genome)

                            ct_mtx_raw[feature_motif_idx, self.bc_mapping[bc], rel_pos] += 1
                            ct_mtx_corr[feature_motif_idx, self.bc_mapping[bc], rel_pos] += (1 / self.bias_correction_factors[kmer])


        self.q.put((ct_mtx_raw, ct_mtx_corr))
        self.sem.release()



def normalize_by_flanks(mtx, flank_size=50):

    mtx_norm = np.zeros(mtx.shape)

    for i in range(mtx.shape[0]):

        flank_avg = np.mean(np.concatenate([mtx[i,:,:flank_size], mtx[i,:,-flank_size:]], axis=1), axis=1)
        mtx_norm[i] = mtx[i] / flank_avg[:, np.newaxis]
        
    return mtx_norm



def smooth_mtx(mtx, window_size=5):

    mtx_smooth = np.zeros(mtx.shape)

    for i in range(mtx.shape[0]):

        for j in range(mtx.shape[2]):

            start = max(j - window_size // 2, 0)
            end = min(j + window_size // 2 + 1, mtx.shape[2])

            mtx_smooth[i, :, j] = np.mean(mtx[i, :, start:end], axis=1)

    return mtx_smooth



def get_scores(adata, mtx):

    center = adata.uns["footprinting"][mtx][0].shape[1] // 2

    tmp = []
    for i in range(len(adata.uns["footprinting"][mtx])):

        motif = adata.uns["motif_search"]["tf_motifs"][i]

        flank_l = adata.uns["footprinting"][mtx][i][:, center - 40 : center - int(np.floor(motif.length / 2))]
        flank_r = adata.uns["footprinting"][mtx][i][:, center + int(np.ceil(motif.length / 2)) : center + 40]
        flank = np.concatenate((flank_l, flank_r), axis=1)

        footprint = adata.uns["footprinting"][mtx][i][:, center - int(np.floor(motif.length / 2)) : center + int(np.ceil(motif.length / 2))]

        # score = np.mean(flank, axis=1) - np.mean(footprint, axis=1)
        score = np.mean(flank, axis=1) / np.mean(footprint, axis=1)

        tmp.append(score)

    score_mtx = np.column_stack(tmp)

    return score_mtx



def footprinting(adata, groupby, fragments, genome, kmer_size=6, distance_to_center=250, smoothing_window_size=None, background="genome", n_jobs=-1):
    
    # create features for the footprinting from the motif search results
    features = get_features(adata, distance_to_center=distance_to_center)

    # load genome
    genome = SeqIO.to_dict(SeqIO.parse(genome, format="fasta"))

    # compute background kmer frequencies
    if background == "peaks":  # peaks as backround
        background_kmer_counts = peak_bg_kmers(adata, kmer_size=kmer_size, genome=genome)
    elif background == "genome":  # genome as background
        background_kmer_counts = genome_bg_kmers(genome, kmer_size=kmer_size, n_jobs=n_jobs)

    background_total = np.sum([count for count in background_kmer_counts.values()])
    kmer_distr_background = {kmer: kmer_count/background_total for kmer, kmer_count in background_kmer_counts.items()}
    print(f"Total background kmers: {background_total}")

    # count number of fragments to prepare for multiprocessing
    n_lines = count_lines(fragments)

    # compute cutsite kmer frequencies
    cutsite_kmer_counts = cutsite_kmers(fragments, genome, adata, n_lines, kmer_size=kmer_size, n_jobs=n_jobs)
    cutsites_total = np.sum([count for count in cutsite_kmer_counts.values()])
    kmer_distr_cutsites = {kmer: kmer_count/cutsites_total for kmer, kmer_count in cutsite_kmer_counts.items()}
    print(f"Total Cutsite kmers: {cutsites_total}")

    # compute bias correction factors
    bias_correction_factors = {kmer: kmer_distr_cutsites[kmer]/kmer_distr_background[kmer] if kmer in kmer_distr_cutsites else 1.0 for kmer in kmer_distr_background.keys()}

    # compute footprinting profiles
    pos_to_feat = build_tree_parallel(features, n_jobs=-1)
    # print(f"{get_total_memory_usage(pos_to_feat) / 1e6} MB")
    tn5_activity(adata, fragments, features, groupby, pos_to_feat, genome, kmer_size, bias_correction_factors, n_lines, region_size=(distance_to_center*2)+1, n_jobs=n_jobs)
    # tn5_activity_parallel(adata, fragments, features, motifs, bias_correction_factors, groupby, genome, aggregate=True, k=k, n_jobs=n_jobs)

    # Normalize (n_obs per cluster)
    # norm_factors = np.array(adata.obs["anno"].value_counts(sort=False))[:, np.newaxis]
    # adata.uns["footprinting"]["mtx_norm"] = [mtx / norm_factors for mtx in adata.uns["footprinting"]["mtx_raw"]]
    # adata.uns["footprinting"]["mtx_corr_norm"] = [mtx / norm_factors for mtx in adata.uns["footprinting"]["mtx_corr"]]

    # flank normalization
    adata.uns["footprinting"]["mtx_norm"] = normalize_by_flanks(adata.uns["footprinting"]["mtx_raw"])
    adata.uns["footprinting"]["mtx_corr_norm"] = normalize_by_flanks(adata.uns["footprinting"]["mtx_corr"])

    # smoothing
    if smoothing_window_size is not None:
        adata.uns["footprinting"]["mtx_norm_smooth"] = smooth_mtx(adata.uns["footprinting"]["mtx_norm"], window_size=smoothing_window_size)
        adata.uns["footprinting"]["mtx_corr_norm_smooth"] = smooth_mtx(adata.uns["footprinting"]["mtx_corr_norm"], window_size=smoothing_window_size)

    # compute scores
    adata.uns["footprinting"]["score_raw"] = get_scores(adata, mtx="mtx_raw")
    adata.uns["footprinting"]["score_norm"] = get_scores(adata, mtx="mtx_norm")
    adata.uns["footprinting"]["score_corr"] = get_scores(adata, mtx="mtx_corr")
    adata.uns["footprinting"]["score_corr_norm"] = get_scores(adata, mtx="mtx_corr_norm")
    if smoothing_window_size is not None:
        adata.uns["footprinting"]["score_norm_smooth"] = get_scores(adata, mtx="mtx_norm_smooth")
        adata.uns["footprinting"]["score_corr_norm_smooth"] = get_scores(adata, mtx="mtx_corr_norm_smooth")



def plot_footprints(adata, mtx="corr_norm", show_score=True, score_threshold=1.3, figsize=None, save=None):

    if mtx in ["raw", "norm", "corr", "corr_norm", "norm_smooth", "corr_norm_smooth"]:
        score = f"score_{mtx}"
        mtx = f"mtx_{mtx}"
    else:
        raise ValueError(f"Unknown mtx '{mtx}'")

    nrows = adata.uns["footprinting"][mtx][0].shape[0]
    ncols = len(adata.uns["footprinting"][mtx])

    if figsize is None:
        figsize = (5 * ncols, 5 * nrows)

    fig, axs = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols, sharey="col")

    x_len = adata.uns["footprinting"][mtx][0].shape[1]
    x = [int(val - ((x_len - 1) / 2)) for val in range(x_len)]


    for i in range(ncols):

        motif_size = adata.uns["motif_search"]["tf_motifs"][i].length
        tf_name = adata.uns["motif_search"]["tf_motifs"][i].name
        mtx_id = adata.uns["motif_search"]["tf_motifs"][i].matrix_id

        min_val = np.inf
        max_val = -np.inf

        for j in range(nrows):

            axs[j, i].plot(x, adata.uns["footprinting"][mtx][i][j], linewidth=1.4)

            if show_score:
                current_score = adata.uns["footprinting"][score][j][i]
                facecolor = "white" if current_score < score_threshold else "green"
                axs[j, i].text(0.05, 0.95, f"{current_score:.2f}", transform=axs[j, i].transAxes, fontsize=12, verticalalignment="top", bbox=dict(boxstyle="round", facecolor=facecolor, alpha=0.75))

            current_min_val = np.min(adata.uns["footprinting"][mtx][i][j])
            if current_min_val < min_val:
                min_val = current_min_val

            current_max_val = np.max(adata.uns["footprinting"][mtx][i][j])
            if current_max_val > max_val:
                max_val = current_max_val

            if i == 0:
                axs[j, i].set_ylabel(adata.uns["footprinting"]["group_mapping"][j], rotation=0, labelpad=50, fontsize="large", fontweight="bold")
                axs[j, i].yaxis.label.set_horizontalalignment("right")

        for j in range(nrows):

            x_motif = [int(i) for i in range(int(-np.floor(motif_size / 2)), int(np.ceil(motif_size / 2)))]

            y_range = max_val - min_val
            y_val_motif = np.max(min_val - (y_range * 0.05), 0)
            y_motif = [y_val_motif for _ in range(len(x_motif))]

            axs[j, i].plot(x_motif, y_motif, linewidth=2, color="tab:red")

        axs[0, i].set_title(f"{tf_name} ({mtx_id}) ({motif_size} bp)", pad=50, fontsize="large", fontweight="bold")
        axs[-1, i].set_xlabel("Distance to motif center [bp]")

    plt.tight_layout()

    if save is not None:
        plt.savefig(save)



def plot_footprints_compact(adata, mtx="corr_norm", max_cols=5, figsize=None, save=None):

    if mtx in ["raw", "norm", "corr", "corr_norm", "norm_smooth", "corr_norm_smooth"]:
        mtx = f"mtx_{mtx}"
    else:
        raise ValueError(f"Unknown mtx '{mtx}'")

    ngroups = adata.uns["footprinting"][mtx][0].shape[0]

    nsubplots = len(adata.uns["footprinting"][mtx])

    ncols = max_cols if nsubplots > max_cols else nsubplots
    nrows = np.ceil(nsubplots / max_cols).astype(int)

    if figsize is None:
        figsize = (5 * ncols, 5 * nrows)

    fig, axs = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols)

    plt.subplots_adjust(hspace=0.4, wspace=0.2)

    axs = axs.flatten()

    for ax in axs[nsubplots+1:]:
        ax.remove()

    axs = axs[:nsubplots+1]

    x_len = adata.uns["footprinting"][mtx][0].shape[1]
    x = [int(val - ((x_len - 1) / 2)) for val in range(x_len)]

    for i, ax in enumerate(axs[:-1]):

        motif_size = adata.uns["motif_search"]["tf_motifs"][i].length
        tf_name = adata.uns["motif_search"]["tf_motifs"][i].name
        mtx_id = adata.uns["motif_search"]["tf_motifs"][i].matrix_id

        min_val = np.inf
        max_val = -np.inf

        for j in range(ngroups):

            if i == len(axs) - 2:
                ax.plot(x, adata.uns["footprinting"][mtx][i][j], linewidth=1.4, label=adata.uns["footprinting"]["group_mapping"][j])
            else:
                ax.plot(x, adata.uns["footprinting"][mtx][i][j], linewidth=1.4)

            current_min_val = np.min(adata.uns["footprinting"][mtx][i][j])
            if current_min_val < min_val:
                min_val = current_min_val

            current_max_val = np.max(adata.uns["footprinting"][mtx][i][j])
            if current_max_val > max_val:
                max_val = current_max_val

        x_motif = [int(i) for i in range(int(-np.floor(motif_size / 2)), int(np.ceil(motif_size / 2)))]

        y_range = max_val - min_val
        y_val_motif = np.max(min_val - (y_range * 0.05), 0)
        y_motif = [y_val_motif for _ in range(len(x_motif))]

        ax.plot(x_motif, y_motif, linewidth=2, color="tab:red")

        ax.set_title(f"{tf_name} ({mtx_id}) ({motif_size} bp)", fontsize="large")
        ax.set_xlabel("Distance to motif center [bp]")

    
    ax_legend = axs[-1]
    ax_legend.axis("off")

    handles, labels = ax.get_legend_handles_labels()
    ax_legend.legend(handles, labels, loc="center", frameon=False)

    if save is not None:
        plt.savefig(save)