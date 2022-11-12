import kneed
import numpy as np
import matplotlib.pyplot as plt


def find_elbow(adata, use_log=False, save=None):
    y = adata.uns["pca"]["variance_ratio"]
    y_log = np.log10(y)

    x = [x for x in range(1, len(y)+1)]

    if not use_log:
        kneedle = kneed.KneeLocator(x, y, curve="convex", direction="decreasing")
    else:
        kneedle = kneed.KneeLocator(x, y_log, curve="convex", direction="decreasing")

    elbow_point = kneedle.elbow

    fig, axs = plt.subplots(figsize=(20,8), nrows=1, ncols=2)
    axs = axs.flatten()

    lwidth = 1.25

    axs[0].plot(x, y, c="black", linewidth=lwidth)
    axs[0].plot(x[elbow_point-2:elbow_point+1], y[elbow_point-2:elbow_point+1], c="tab:green", linewidth=lwidth)
    axs[0].plot(elbow_point, y[elbow_point-1], marker="o", c="tab:green", markersize=12)
    axs[0].plot(x, y, marker=".", c="black", markersize=6, linewidth=0)

    axs[0].annotate(
        "Elbow point = {}".format(elbow_point),
        xy=(elbow_point, y[elbow_point-1]),
        xytext=(15, 50), textcoords="offset points",
        arrowprops=dict(arrowstyle="-|>", mutation_scale=20, linewidth=2, color="black"),
        ha="left", va="center",
        fontsize=14, fontweight="bold"
    )

    axs[0].set_title("Raw")

    axs[1].plot(x, y_log, c="black", linewidth=lwidth)
    axs[1].plot(x[elbow_point-2:elbow_point+1], y_log[elbow_point-2:elbow_point+1], c="tab:green", linewidth=lwidth)
    axs[1].plot(elbow_point, y_log[elbow_point-1], marker="o", c="tab:green", markersize=12)
    axs[1].plot(x, y_log, marker=".", c="black", markersize=6, linewidth=0)

    axs[1].annotate(
        "Elbow point = {}".format(elbow_point),
        xy=(elbow_point, y_log[elbow_point-1]),
        xytext=(15, 50), textcoords="offset points",
        arrowprops=dict(arrowstyle="-|>", mutation_scale=20, linewidth=2, color="black"),
        ha="left", va="center",
        fontsize=14, fontweight="bold"
    )

    axs[1].set_title("Log Transformed")

    if not save:
        plt.show()

    else:
        if isinstance(save, str):
            filename = save
        else:
            filename = "kneedle.png"

        plt.savefig(filename, dpi=300)

    return elbow_point