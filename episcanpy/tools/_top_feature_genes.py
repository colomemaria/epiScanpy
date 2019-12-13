import warnings
from warnings import warn

def var_features_to_genes(adata, gtf_file, extension=5000):
    """
    Once you called the most variable features.
    You can identify genes neighboring these features of interest.

    """
    
    # extract_top_markers
    #print(adata.uns['rank_genes_groups'].keys())
    windows = [list(w) for w in adata.uns['rank_genes_groups']['names'].tolist()]
    windows_all = []
    for w in windows:
        windows_all += w

    # load the gtf file
    with open(gtf_file) as f:
        gtf_raw = f.readlines()
    gtf = []
    for line in gtf_raw:
        if line[0] != '#':
            gtf.append(line[:-2].split('\t'))
    del gtf_raw
    
    gtf_dic = {}
    for line in gtf:
        if line[0] not in gtf_dic.keys():
            gtf_dic[line[0]] = [line]
        else:
            gtf_dic[line[0]].append(line)
    del gtf

    markers = []
    for w in windows_all:
        curr_m = []
        w2 = w.split('_')
        chrom = w2[0][3:]
        w2 = [int(x) for x in w2[1:]]
        for gene in gtf_dic[chrom]:
            start = int(gene[3])-extension
            end = int(gene[4])+extension
            if (w2[0] < start < w2[1]) or (w2[0] < end < w2[1]):
                gene_name = gene[-1]
                for n in gene[-1].split(';'):
                    if 'gene_name' in n:
                        gene_name = n
                curr_m.append([w, gene_name])
        if curr_m != []:
            markers += curr_m
    markers = [list(x) for x in set(tuple(x) for x in markers)]
    markers_dict={}
    for n in markers:
        markers_dict[n[1].split(' "')[1][:-1]] = n[0]
    return(markers_dict)


def top_feature_genes(adata, gtf_file, extension=5000):
    """
    Deprecated - Please use epi.tl.var_features_to_genes instead.
    Once you called the most variable features.
    You can identify genes neighboring these features of interest.

    """
    warn.warn('Deprecated - Please use epi.tl.var_features_to_genes instead.')
    var_features_to_genes(adata=adata, gtf_file=gtf_file, extension=extension)