import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.axes as pltax
import numpy as np
import anndata as ad

import warnings
from warnings import warn

from scipy.sparse import issparse
from scipy.stats.stats import pearsonr, spearmanr

def cal_var(adata, show=True, color=['b', 'r'], save=None):
    """
    Show distribution plots of cells sharing features and variability score.

    Parameters
    ----------

    adata

    show

    color

    save

    """

    if issparse(adata.X):
        adata.var['n_cells'] = adata.X.sum(axis=0).tolist()[0] ## double check this
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = [1-abs(n-0.5) for n in adata.var['prop_shared_cells']]
    else:
        adata.var['n_cells'] = adata.X.sum(axis=0)
        adata.var['prop_shared_cells'] = adata.var['n_cells']/len(adata.obs_names.tolist())
        adata.var['variability_score'] = [1-abs(n-0.5) for n in adata.var['prop_shared_cells']]

    if save!= None:
        fig, axs = plt.subplots(ncols=2)
        plt.subplots_adjust(wspace=0.7)

        ax0 = sns.distplot(adata.var['prop_shared_cells'], bins=40, ax=axs[0], color=color[0])
        ax0.set(xlabel='cells sharing a feature', ylabel='density')
        ax1 = sns.distplot(adata.var['variability_score'], bins=40, ax=axs[1], color=color[1])
        ax1.set(xlabel='variability score', ylabel='density')
        plt.savefig(save, bbox_inches="tight")
        show = False

    if show: # plotting
        fig, axs = plt.subplots(ncols=2)
        plt.subplots_adjust(wspace=0.7)

        ax0 = sns.distplot(adata.var['prop_shared_cells'], bins=40, ax=axs[0], color=color[0])
        ax0.set(xlabel='cells sharing a feature', ylabel='density')
        ax1 = sns.distplot(adata.var['variability_score'], bins=40, ax=axs[1], color=color[1])
        ax1.set(xlabel='variability score', ylabel='density')
        plt.show()

def select_var_feature(adata, min_score=0.5, nb_features=None, show=True, copy=False, save=None):
    """
    This function computes a variability score to rank the most variable features across all cells.
    Then it selects the most variable features according to either a specified number of features (nb_features) or a maximum variance score (max_score).

    Parameters
    ----------

    adata: adata object

    min_score: minimum threshold variability score to retain features,
    where 1 is the score of the most variable features and 0.5 is the score of the least variable features.

    nb_features: default value is None, if specify it will select a the top most variable features.
    if the nb_features is larger than the total number of feature, it filters based on the max_score argument

    show: default value True, it will plot the distribution of var.

    copy: return a new adata object if copy == True.

    Returns
    -------
    Depending on ``copy``, returns a new AnnData object or overwrite the input


    """
    if copy:
        inplace=False
    else:
        inplace=True

    adata = adata.copy() if not inplace else adata

    # calculate variability score
    cal_var(adata, show=show, save=save) # adds variability score for each feature
    # adata.var['variablility_score'] = abs(adata.var['prop_shared_cells']-0.5)
    var_annot = adata.var.sort_values(ascending=False, by ='variability_score')

    # calculate the max score to get a specific number of feature
    if nb_features != None and nb_features < len(adata.var_names):
            min_score = var_annot['variability_score'][nb_features]


    adata_tmp = adata[:,adata.var['variability_score']>=min_score].copy()

    ## return the filtered AnnData objet.
    if not inplace:
        adata_tmp = adata[:,adata.var['variability_score']>=min_score]
        return(adata_tmp)
    else:
        adata._inplace_subset_var(adata.var['variability_score']>=min_score)


def binarize(adata, copy=False):
    """convert the count matrix into a binary matrix.

    Parameters
    ----------

    adata: AnnData object

    copy: return a new adata object if copy == True.

    Returns
    -------
    Depending on ``copy``, returns a new AnnData object or overwrite the input

    """

    if copy:
        adata2 = adata.copy()
        adata2.X[adata2.X != 0] = 1
        return(adata2)
    else:
        adata.X[adata.X != 0] = 1



def coverage_cells(adata,
                   key_added=None,
                   log=False,
                   binary=None,
                   ## threshold
                   threshold=None,
                   bw=0.5,
                   ## plotting
                   bins=50,
                   xlabel=None, ylabel=None, title=None,
                   #figsize=(15,10),
                   #save_dpi=250,
                   color=None, edgecolor=None,
                   ## sving
                   save=None):
    """
    Histogram of the number of open features (in the case of ATAC-seq data) per cell.

    Parameters
    ----------

    adata

    binary
    log : default not . Accepted log2, log10, log1p, log (natural log base e).
    if log is True, log10 used

    bins
    threshold

    key_added

    xlabel
    ylabel
    title
    color
    edgecolor
    save

    """
    if key_added == None:
        key_added='nb_features'

    # calculate the number of features per cell
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        sum_peaks = np.sum(adata.X, axis=1).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        sum_peaks = np.sum(tmp_array.X, axis=1).tolist()

    if issparse(adata.X):
        sum_peaks = [element[0] for element in sum_peaks]

    adata.obs[key_added] = sum_peaks

    #fig = plt.figure(figsize=figsize)
    # plotting settings
    if xlabel ==None:
        plt.xlabel('number of open features per cell ')
    else:
        plt.xlabel(xlabel)

    if ylabel ==None:
        plt.ylabel('number of cells')
    else:
        plt.ylabel(ylabel)

    if title !=None:
        plt.title(title)

    if color == None:
        color='c'
    if edgecolor == None:
        edgecolor='k'

    bw_param = bw
    sns.set_style('whitegrid')

    ### Potential log scale
    if log!=False:
        if 0 in sum_peaks:
            warnings.warn("""Some cells do not contain any open feature. Use epi.pp.filter_cells(adata, max_features=1) to remove these cells.""")

        if log=='log2':
            plt.xlabel('number of open features per cell (log2)')
            fig = plt.hist(np.log2(sum_peaks), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log2(threshold), color='r', linestyle='--')
        elif log=='log1p':
            plt.xlabel('number of open features per cell (log1p)')
            fig = plt.hist(np.log1p(sum_peaks), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log1p(threshold), color='r', linestyle='--')
        elif log=='log':
            plt.xlabel('number of open features per cell (log natural base e)')
            fig = plt.hist(np.log(sum_peaks), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log(threshold), color='r', linestyle='--')
        else:
            plt.xlabel('number of open features per cell (log10)')
            fig = plt.hist(np.log10(sum_peaks), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log10(threshold), color='r', linestyle='--')

    else:
        fig = plt.hist(sum_peaks, bins, color=color, edgecolor=edgecolor)
        if threshold != None:
            plt.axvline(x=threshold, color='r', linestyle='--')


    if save!= None:
        #fig.savefig(save, dpi=save_dpi)
        plt.savefig(save, bbox_inches="tight")
    plt.show()
    adata.obs[key_added] = sum_peaks

def coverage_features(adata,
                        binary=None,
                        log=False,
                        key_added=None,
                        threshold=None,
                        bw=0.5,
                        bins=50,
                        xlabel=None, ylabel=None, title=None,
                        #figsize=(15,10),
                        #save_dpi=250,
                        color=None, edgecolor=None,
                        save=None):

    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.

    Parameters
    ----------
    adata
    threshold
    log
    binary
    key_added
    bw
    xlabel
    title
    color
    edgecolor

    """

    if key_added == None:
        key_added='commonness'

    # calculate for each feature in how many cell it is open
    if binary==None:
        warnings.warn("""The argument binary was not specified. To reduce computing time, you can specify if the matrix is already binary""")
    if binary:
        common = np.sum(adata.X, axis=0).tolist()
    else:
        tmp_array = binarize(adata, copy=True)
        common = np.sum(tmp_array.X, axis=0).tolist()

    if issparse(adata.X):
        common = common[0]
    adata.var[key_added] = common

    #fig = plt.figure(figsize=figsize)
    ### plotting settings
    if xlabel ==None:
        plt.xlabel('number of cells sharing a feature')
    else:
        plt.xlabel(xlabel)
    if ylabel ==None:
        plt.ylabel('number of features')
    else:
        plt.ylabel(ylabel)
    if title !=None:
        plt.title(title)

    bw_param = bw
    sns.set_style('whitegrid')

    if color == None:
        color='c'
    if edgecolor == None:
        edgecolor='k'


    ### Potential log scale
    if log !=False:
        if 0 in common:
            warnings.warn("""Some features are not open in any cell. Use epi.pp.filter_features(adata, max_cells=1) to remove these features.""")

        # different potential bases for log transfo
        if log=='log2':
            plt.xlabel('number of cells sharing a feature (log2)')
            fig = plt.hist(np.log2(common), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log2(threshold), color='r', linestyle='--')
        elif log=='log1p':
            plt.xlabel('number of cells sharing a feature (log1p)')
            fig = plt.hist(np.log1p(common), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log1p(threshold), color='r', linestyle='--')
        elif log=='log':
            plt.xlabel('number of cells sharing a feature (log natural base e)')
            fig = plt.hist(np.log(common), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log(threshold), color='r', linestyle='--')
        else:
            plt.xlabel('number of cells sharing a feature (log10)')
            fig = plt.hist(np.log10(common), bins, color=color, edgecolor=edgecolor)
            if threshold != None:
                plt.axvline(x=np.log10(threshold), color='r', linestyle='--')

    else:
        fig = plt.hist(common, bins=int(80), color=color, edgecolor=edgecolor)
        if threshold != None:
                plt.axvline(x=threshold, color='r', linestyle='--')


    if save!= None:
        #fig.savefig(save, dpi=save_dpi)
        plt.savefig(save, bbox_inches="tight")
    plt.show()

    adata.var[key_added] = common





def correlation_pc(adata,
                   variable,
                   pc=1,
                   obs=True,
                   method='pearson',
                   title=None,
                   xlabel=None,
                   ylabel=None,
                   #figsize=(15,10),
                   #save_dpi=250,
                   show=True,
                   save=None,
                  ):
    """
    Correlation between a given PC and a covariate.
    If show == True, plot a scatter plot.
    If obs == True, consider a obs covariate. Else, take a variable covariate.
    Available methods for correlation: 'pearson' or 'spearman'

    Parameters
    ----------

    adata : input adata

    variable : covariate either saved in obs or var

    obs : Boolean that specify if the covariate is stored in obs or in var.
    In later version, this parameter will not need to be specified

    method : 'pearson' or 'spearman' available - scipy.stats implementation

    title : optional title to the plot

    xlabel : optional xlabel

    ylabel : optional ylabel

    show: Print the correlation coefficient and p-value. Additionaly, plot a scatter plot
    of the PC coordinate and the covariate value for every cell or feature in your matrix

    save: if specified the name of the picture, save as png

    Return
    ------

    correlation coefficient and p-value

    """
    if pc-1>len(adata.varm['PCs'][0]):
        warnings.warn("".join(["""You requested a PC that is not currently available.
                            Please run epi.pp.pca(adata, n_comps=""", int(pc+1), ') ']))


    if obs:
        x = adata.obsm['X_pca'][:,pc-1]
        y = adata.obs[variable]
    else:
        x = adata.varm['PCs'][:,pc-1]
        y = adata.var[variable]

    if method=='pearson':
        correlation = pearsonr(x, y)
    elif method=='spearman':
        correlation = spearmanr(x, y)

    #fig = plt.figure(figsize=figsize)
    if show:
        ## Add legends and title
        if xlabel:
            plt.xlabel(xlabel)
        else:
            plt.xlabel(" ".join(['X_PC', str(pc)]))

        if ylabel:
            plt.ylabel(ylabel)
        else:
            plt.ylabel(variable)

        if title:
            plt.title(title)

        plt.scatter(x, y)
        print("correlation: " +str(correlation[0]))
        print("pval: " +str(correlation[1]))

    adata.uns['correlation_pc'] = correlation

    if save!= None:
        #fig.savefig(save, dpi=save_dpi)
        plt.savefig(save, bbox_inches="tight")
    plt.show()


def commonness_features(adata,
                        binary=None,
                        log=False,
                        key_added=None,

                        threshold=None,
                        bw=0.5,
                        bins=50,
                        xlabel=None, ylabel=None, title=None,
                        color=None, edgecolor=None,

                        save=None):

    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.

    Parameters
    ----------
    adata
    threshold
    log
    binary
    key_added
    bw
    xlabel
    title
    color
    edgecolor

    """
    warnings.warn(""" Function deprecated. Use epi.pp.coverage_features instead.
        Or use epi.pp.density_features.
    """)

    coverage_features(adata,
                        binary=binary,
                        log=log,
                        key_added=key_added,
                        threshold=threshold,
                        bw=bw,
                        bins=bins,
                        xlabel=xlabel, ylabel=ylabel, title=title,
                        color=color, edgecolor=edgecolor,
                        save=save)


def density_features(adata,
    threshold=None,
    hist=False,
    bins=80,
    bw=0.5,
    key_added=None,
    xlabel=None,
    title=None,
    figsize=(15,10),
    save_dpi=250,
    save=None):
    """
    Display how often a feature is measured as open (for ATAC-seq).
    Distribution of the feature commoness in cells.
    """
    if key_added == None:
        key_added='commonness'

    common = np.sum(adata.X, axis=0).tolist()
    if len(common) == 1:
        common = [item for sublist in common for item in sublist]
    adata.var[key_added] = common


    fig = plt.figure(figsize=figsize)
    if threshold != None:
        plt.axvline(x=threshold, color='r')

    bw_param = bw
    sns.set_style('whitegrid')
    sns.distplot(common, hist=hist, kde=True,
                    bins=bins, color = 'darkblue',
                    hist_kws={'edgecolor':'black'},
                    kde_kws={'linewidth': 1})
    if xlabel ==None:
        plt.xlabel('cells sharing a feature')
    else:
        plt.xlabel(xlabel)

    plt.ylabel('density')

    if title !=None:
        plt.title(title)


    if save!= None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()

    adata.var[key_added] = common




def variability_features(adata, min_score=None, nb_features=None, show=True,
                         title=None,
                         log=None,
                         xlabel='ranked features',
                         ylabel='variability score',
                         #figsize=(15,10),
                         #save_dpi=250,
                         save=None):
    """
    This function computes a variability score to rank the most variable features across all cells.
    Then it selects the most variable features according to either a specified number of features (nb_features) or a maximum variance score (max_score).

    Parameters
    ----------

    adata: adata object

    min_score:  default value is None, if specify it will plot threshold

    nb_features: default value is None, if specify it will plot threshold

    log: default false. You can specify 'log2' or 'log10' to put the ranked features on a log scale

    show: default value True

    save: if specified, save the figure.

    """
    #fig = plt.figure(figsize=figsize)
    # calculate variability score
    cal_var(adata, show=False) # adds variability score for each feature
    # sort features per variability score
    var_annot = adata.var.sort_values(ascending=False, by ='variability_score')
    all_feature = [n for n in range(1, len(var_annot)+1)]
    # add log scale
    if log != None:
        if log=='log2':
            all_feature = np.log2(all_feature)
            xlabel = xlabel + ' (log2)'
            if nb_features != None:
                nb_features = np.log2(nb_features)
        elif log=='log10':
            all_feature = np.log10(all_feature)
            xlabel = xlabel + ' (log10)'
            if nb_features != None:
                nb_features = np.log10(nb_features)
        else:
            warn("invalid log option. Accepted log scales are 'log2' and 'log10'.")

    # add min_score threshold
    if min_score != None:
        plt.axhline(min_score, color='r', linestyle='--')

    # add nb_features threshold
    if nb_features != None:
        plt.axvline(nb_features, color='r', linestyle='--')

    plt.plot(all_feature, var_annot.variability_score)
    # title and axis
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if save != None:
        #fig.savefig(save, dpi=save_dpi)
        plt.savefig(save, bbox_inches="tight")
    plt.show()
    #return(var_annot)


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

    Args:
        adata: AnnData
        gtf: path to GTF file
        fragments: path to fragments file
        method: method to use for computing TSS enrichment
        score: value that is used as TSS enrichment score for individual observations
        distance_to_tss: distance to TSS
        bp_per_flank: number of base pair per flank to use for normalization
        n_jobs: number of jobs to use for computation. -1 means using all processors

    Returns:
        None
    """

    if method not in {"standard", "complete_coverage"}:
        raise ValueError("value for method argument must be 'standard' or 'complete_coverage'")

    if score not in {"avg_score_of_center_region", "score_at_zero"}:
        raise ValueError("value for score argument must be 'avg_score_of_center_region' or 'score_at_zero'")

    import multiprocessing
    import numpy as np
    import pandas as pd
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


def tss_enrichment_plot(adata,
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

    import matplotlib.pyplot as plt
    import numpy as np

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
        default_filename = "tss_enrichment_score"

        if isinstance(save, str):

            if save[-4:] in {".pdf", ".png", ".svg"}:
                filename = default_filename + save

            else:
                filename = default_filename + save + ".png"

        else:
            filename = default_filename + ".png"

        plt.savefig(filename)


def tss_enrichment_score_plot(adata, figsize=(4, 6), save=None):
    """
    Plots a violin plot of the individual TSS enrichment scores

    Args:
        adata: AnnData
        figsize: size of the figure
        save: if True or str, save the figure. str is appended to the default filename. infer filetype if ending on {'.pdf', '.png', '.svg'}

    Returns:
        None
    """

    import seaborn as sns
    import matplotlib.pyplot as plt

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
        default_filename = "tss_enrichment_score"

        if isinstance(save, str):

            if save[-4:] in {".pdf", ".png", ".svg"}:
                filename = default_filename + save

            else:
                filename = default_filename + save + ".png"

        else:
            filename = default_filename + ".png"

        plt.savefig(filename)
