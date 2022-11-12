import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def violin(adata, key, min_threshold=None, max_threshold=None, show_log=True, show_mean=True, show_median=True, print_statistics=True, save=None):

    figsize = (6, 6) if not show_log else (12, 6)
    ncols = 1 if not show_log else 2

    fig, axs = plt.subplots(figsize=figsize, nrows=1, ncols=ncols, squeeze=False)
    axs = axs.flatten()

    df = adata.obs if key in adata.obs else adata.var

    sns.violinplot(y=df[key], inner=None, cut=0, ax=axs[0])

    if show_mean:
        axs[0].axhline(y=df[key].mean(), color="white", linestyle="--", linewidth=2, alpha=0.75)

    if show_median:
        axs[0].axhline(y=df[key].median(), color="white", linestyle="-", linewidth=2, alpha=0.75)

    if max_threshold and min_threshold:
        sns.stripplot(y=df[key][np.logical_and(df[key] <= max_threshold, df[key] >= min_threshold)], size=1, jitter=0.4, color="black", ax=axs[0])
    elif max_threshold:
        sns.stripplot(y=df[key][df[key] <= max_threshold], size=1, jitter=0.4, color="black", ax=axs[0])
    elif min_threshold:
        sns.stripplot(y=df[key][df[key] >= min_threshold], size=1, jitter=0.4, color="black", ax=axs[0])
    else:
        sns.stripplot(y=df[key], size=1, jitter=0.4, color="black", ax=axs[0])

    if max_threshold:
        sns.stripplot(y=df[key][df[key] > max_threshold], size=1, jitter=0.4, color="red", ax=axs[0])
        axs[0].axhline(y=max_threshold, color="red", linestyle="--", linewidth=1, alpha=0.75)
    if min_threshold:
        sns.stripplot(y=df[key][df[key] < min_threshold], size=1, jitter=0.4, color="red", ax=axs[0])
        axs[0].axhline(y=min_threshold, color="red", linestyle="--", linewidth=1, alpha=0.75)


    if show_log:
        sns.violinplot(y=[np.log10(val) if val != 0 else 0 for val in df[key].values], inner=None, cut=0, ax=axs[1])

        if show_mean:
            axs[1].axhline(y=np.log10(df[key].mean()), color="white", linestyle="--", linewidth=2, alpha=0.75)

        if show_median:
            axs[1].axhline(y=np.log10(df[key].median()), color="white", linestyle="-", linewidth=2, alpha=0.75)

        if max_threshold and min_threshold:
            sns.stripplot(y=np.log10(df[key][np.logical_and(df[key] <= max_threshold, df[key] >= min_threshold)]), size=1, jitter=0.4, color="black", ax=axs[1])
        elif max_threshold:
            sns.stripplot(y=np.log10(df[key][df[key] <= max_threshold]), size=1, jitter=0.4, color="black", ax=axs[1])
        elif min_threshold:
            sns.stripplot(y=np.log10(df[key][df[key] >= min_threshold]), size=1, jitter=0.4, color="black", ax=axs[1])
        else:
            sns.stripplot(y=np.log10(df[key]), size=1, jitter=0.4, color="black", ax=axs[1])

        if max_threshold:
            sns.stripplot(y=np.log10(df[key][df[key] > max_threshold]), size=1, jitter=0.4, color="red", ax=axs[1])
            axs[1].axhline(y=np.log10(max_threshold), color="red", linestyle="--", linewidth=1, alpha=0.75)
        if min_threshold:
            sns.stripplot(y=np.log10(df[key][df[key] < min_threshold]), size=1, jitter=0.4, color="red", ax=axs[1])
            axs[1].axhline(y=np.log10(min_threshold), color="red", linestyle="--", linewidth=1, alpha=0.75)

        axs[1].set_ylabel("log10({})".format(key))

    plt.tight_layout()

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "qc_violin_{}.png".format(key)

        plt.savefig(filename, dpi=300)

    if print_statistics:
        print("Max:\t{}".format(df[key].max()))
        print("Median:\t{}".format(df[key].median()))
        print("Mean:\t{}".format(df[key].mean()))
        print("Min:\t{}".format(df[key].min()))
