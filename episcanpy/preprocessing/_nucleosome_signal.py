import numpy as np
import gzip
import matplotlib.pyplot as plt
import seaborn as sns

def nucleosome_signal(adata, fragments, n=10000):

    n = adata.n_obs * n

    valid_bcs = set(adata.obs.index)

    comment = "#"
    check_for_comments = True
    use_strip = None

    fragment_info = {bc: [0, 0, []] for bc in valid_bcs}

    if fragments.endswith(".gz"):
        fh = gzip.open(fragments, mode="rt")
    else:
        fh = open(fragments, mode="r")

    c = 0

    for line in fh:

        if check_for_comments:
            if line.startswith(comment):
                continue
            else:
                check_for_comments = False

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

        bc = line_split[3]
        if bc not in valid_bcs:
            continue

        chrom = line_split[0]
        start = int(line_split[1])
        stop = int(line_split[2])

        length = stop - start
        if length < 147:
            fragment_info[bc][0] += 1
        elif length >= 147 and length < 294:
            fragment_info[bc][1] += 1

        fragment_info[bc][2].append(length)

        c += 1
        if c == n:
            break

    nucleosome_signal = dict()

    for bc, (nuc_free, mono_nuc, lns) in fragment_info.items():
        tmp = (
            mono_nuc / nuc_free if nuc_free != 0 else np.nan,
            lns
        )
        nucleosome_signal[bc] = tmp

    adata.obs["nucleosome_signal"] = [nucleosome_signal[bc][0] for bc in adata.obs.index]
    adata.uns["fragment_lengths"] = {bc: nucleosome_signal[bc][1] for bc in adata.obs.index}


def fragment_length(adata, n=5000, threshold=4, show_n=True, save=None):

    ncols = 1 if not threshold else 2
    nrows = 1

    fig, axs = plt.subplots(figsize=(ncols*5, nrows*5), nrows=nrows, ncols=ncols, squeeze=False)
    axs = axs.flatten()

    for i, ax in enumerate(axs):

        if not threshold:
            vals = [[length for length in adata.uns["fragment_lengths"][bc] if length <= 800][:n] for bc in adata.obs.index]
            n_obs = adata.n_obs
            color = "tab:blue"
            vals = [y for x in vals for y in x[:n]]

        else:
            if i == 0:
                vals = [[length for length in adata.uns["fragment_lengths"][bc] if length <= 800][:n] for bc in adata[adata.obs.nucleosome_signal <= threshold].obs.index]
                n_obs = adata[adata.obs.nucleosome_signal <= threshold].n_obs
                color = "tab:blue"
            if i == 1:
                vals = [[length for length in adata.uns["fragment_lengths"][bc] if length <= 800][:n] for bc in adata[adata.obs.nucleosome_signal > threshold].obs.index]
                n_obs = adata[adata.obs.nucleosome_signal > threshold].n_obs
                color = "tab:orange"

            vals = [y for x in vals for y in x[:n]]

        ax.hist(vals,
                bins=200,
                linewidth=0,
                color=color)

        ax.set_xlim((0, 800))

        if show_n:
            ax.annotate(text="n = {}".format(n_obs),
                        xy=(0.8, 0.8),
                        xycoords="axes fraction",
                        xytext=(0, 0),
                        textcoords="offset points",
                        fontsize=12,
                        ha="right",
                        va="bottom")

    sns.despine()
    plt.tight_layout()

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "fragment_length.png"

        plt.savefig(filename, dpi=300)
