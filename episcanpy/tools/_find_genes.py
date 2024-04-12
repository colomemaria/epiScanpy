import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree

def find_genes(adata,
                 gtf_file,
                 key_added='gene_annotation',
                 upstream=5000,
                 downstream=0,
                 feature_type='gene',
                 annotation='HAVANA',
                 raw=False):
    """
    merge values of peaks/windows/features overlapping genebodies + 2kb upstream.
    It is possible to extend the search for closest gene to a given number of bases downstream as well.

    There is commonly 2 set of annotations in a gtf file(HAVANA, ENSEMBL). By default, the function
    will search annotation from HAVANA but other annotation label/source can be specifed.

    It is possible to use other type of features than genes present in a gtf file such as transcripts or CDS.

    """
    ### extracting the genes
    gtf = {}
    with open(gtf_file) as f:
        for line in f:
            if line[0:2] != '##' and '\t'+feature_type+'\t' in line and '\t'+annotation+'\t' in line:
                line = line.rstrip('\n').split('\t')
                if line[6] == '-':
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3])-downstream, int(line[4])+upstream,line[-1].split(';')[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3])-downstream, int(line[4])+upstream,line[-1].split(';')[:-1]])
                else:
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3])-upstream, int(line[4])+downstream,line[-1].split(';')[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3])-upstream, int(line[4])+downstream,line[-1].split(';')[:-1]])

    # extracting the feature coordinates
    raw_adata_features = {}
    feature_index = 0
    for line in adata.var_names.tolist():
        line = line.split('_')
        if line[0] not in raw_adata_features.keys():
            raw_adata_features[line[0]] = [[int(line[1]),int(line[2]), feature_index]]
        else:
            raw_adata_features[line[0]].append([int(line[1]),int(line[2]), feature_index])
        feature_index += 1

    ## find the genes overlaping the features.
    gene_index = []
    for chrom in raw_adata_features.keys():
        if chrom in gtf.keys():
            chrom_index = 0
            previous_features_index = 0
            for feature in raw_adata_features[chrom]:
                gene_name = []
                feature_start = feature[0]
                feature_end = feature[1]
                for gene in gtf[chrom]:
                    if (gene[1]<= feature_start): # the gene is before the feature. we need to test the next gene.
                        continue
                    elif (feature_end <= gene[0]): # the gene is after the feature. we need to test the next feature.
                        break
                    else: # the window is overlapping the gene.
                        for n in gene[-1]:
                            if 'gene_name' in n:
                                gene_name.append(n.lstrip('gene_name "').rstrip('""'))

                if gene_name == []:
                    gene_index.append('intergenic')
                elif len(gene_name)==1:
                    gene_index.append(gene_name[0])
                else:
                    gene_index.append(";".join(list(set(gene_name))))

        else:
            for feature in raw_adata_features[chrom]:
                gene_index.append("unassigned")

    adata.var[key_added] = gene_index



def read_gtf(gtf_file, source=None, gene_type=None):

    names = ["chr", "source", "type", "start", "stop", "score", "strand", "frame", "attribute"]
    annotations = pd.read_csv(gtf_file, sep="\t", header=None, comment="#", names=names, dtype={"chr": str})

    annotations = annotations[annotations.type == "gene"]

    if source:
        annotations = annotations[annotations.source == source]

    annotations["gene_id"] = [attr.replace("gene_id", "").strip().strip("\"") for feature_attr in annotations.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_id")]

    tmp = [attr.replace("gene_name", "").strip().strip("\"") for feature_attr in annotations.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_name ")]
    if not tmp:
        tmp = [attr.replace("gene", "").strip().strip("\"") for feature_attr in annotations.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene ")]
    annotations["gene_name"] = tmp

    tmp = [attr.replace("gene_type", "").strip().strip("\"") for feature_attr in annotations.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_type ")]
    if not tmp:
        tmp = [attr.replace("gene_biotype", "").strip().strip("\"") for feature_attr in annotations.attribute for attr in feature_attr.split(";") if attr.strip().startswith("gene_biotype ")]
    annotations["gene_type"] = tmp

    if gene_type:
        annotations = annotations[[feature in gene_type for feature in annotations.gene_type]]

    annotations.index = annotations.gene_id

    annotations.sort_values(by=["chr", "start", "stop"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)

    annotations = annotations[["gene_name", "gene_id", "gene_type", "chr", "start", "stop", "strand", "source"]]

    return annotations



def add_basal_regulatory_domain(annotations, upstream=5000, downstream=1000):

    annotations["start_brd"] = annotations.apply(lambda row: row.start - upstream if row.strand == "+" else row.start - downstream, axis=1)
    annotations["start_brd"] = annotations["start_brd"].clip(lower=0)
    annotations["stop_brd"] = annotations.apply(lambda row: row.start + downstream if row.strand == "+" else row.start + upstream, axis=1)

    annotations.sort_values(by=["chr", "start_brd", "stop_brd"], key=lambda col: col if col.dtype == np.int64 else col.str.lower(), inplace=True)
    
    return annotations



def add_extended_regulatory_domain(annotations, max_extension=1000000):

    annotations["start_erd"] = annotations.apply(lambda row: row.start - max_extension, axis=1)
    annotations["start_erd"] = annotations["start_erd"].clip(lower=0)
    annotations["stop_erd"] = annotations.apply(lambda row: row.stop + max_extension, axis=1)

    for i in range(1, len(annotations)):

        prev_row = annotations.iloc[i-1]
        current_row = annotations.iloc[i]

        if prev_row.chr != current_row.chr:
            continue

        if prev_row.stop_erd >= current_row.start_brd:
            annotations.at[prev_row.name, "stop_erd"] = np.max([current_row.start_brd - 1, prev_row.stop_brd])
            annotations.at[current_row.name, "start_erd"] = np.min([prev_row.stop_brd + 1, current_row.start_brd])

    return annotations



def add_extended_gene_body(annotations, upstream=5000, downstream=0):

    annotations["start_egb"] = annotations.apply(lambda row: row.start - upstream if row.strand == "+" else row.start - downstream, axis=1)
    annotations["start_egb"] = annotations["start_egb"].clip(lower=0)
    annotations["stop_egb"] = annotations.apply(lambda row: row.stop + downstream if row.strand == "+" else row.stop + upstream, axis=1)

    return annotations



def annotate_peaks(adata, gtf_file, mode="extended_regulatory_domain", source="HAVANA", gene_type=None):

    annotations = read_gtf(gtf_file, source=source, gene_type=gene_type)

    if mode in ["basal_regulatory_domain", "extended_regulatory_domain"]:
        annotations = add_basal_regulatory_domain(annotations, upstream=5000, downstream=1000)

    if mode == "extended_regulatory_domain":
        annotations = add_extended_regulatory_domain(annotations, max_extension=1000000)

    if mode == "extended_gene_body":
        annotations = add_extended_gene_body(annotations, upstream=5000, downstream=0)

    # create a genomic tree
    tree = {chrom: IntervalTree() for chrom in annotations.chr.unique()}
    for annotation in annotations.itertuples():

        if mode == "basal_regulatory_domain":
            interval_start, interval_stop = annotation.start_brd, annotation.stop_brd
        elif mode == "extended_regulatory_domain":
            interval_start, interval_stop = annotation.start_erd, annotation.stop_erd
        elif mode == "extended_gene_body":
            interval_start, interval_stop = annotation.start_egb, annotation.stop_egb
            
        tree[annotation.chr].add(Interval(interval_start, interval_stop, annotation.gene_name))

    # represent feature locations as single base pair positions (left median of the span)
    feature_locations = adata.var.apply(lambda row: (row.chr, np.floor(np.median([row.start, row.stop])).astype(int)), axis=1)
    
    assigned_annotations = []
    for chrom, pos in feature_locations:

        try:
            assigned_annotation = tree[chrom].at(pos)
            assigned_annotation = ";".join([entry.data for entry in assigned_annotation]) if assigned_annotation else "intergenic"
        except KeyError:
            assigned_annotation = "unassigned"
            
        assigned_annotations.append(assigned_annotation)

    adata.var[f"annotation_{mode}"] = assigned_annotations

    return adata