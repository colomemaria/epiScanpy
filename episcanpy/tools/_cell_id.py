import scanpy as sc
import warnings
from collections import Counter

def load_markers(path, marker_list_file):
    """
    Convert list of known cell type markers from literature to a dictionary
    Input list of known marker genes 
    First row is considered the header

    Paramters
    ---------
    marker_list_file: file with the known marker genes (gene name is the 3rd column)
    
    Return
    ------
    cell_type_markers

    """
    marker_list = []
    cell_type_list = []
    with open(path+marker_list_file) as f:
        head = f.readline()
        for line in f:
            line = line.split("\t")
            line[-1] = line[-1][:-1]
            marker_list.append(line[2].upper())
            cell_type_list.append(line[0])
            
    unique_types = list(set(cell_type_list))
    
    # dictionary: list of marker genes per cell type
    cell_type_markers = {}
    for i in range(0,len(set(cell_type_list))):
        cell_type_markers[unique_types[i]] = []
        for j in range(0, len(cell_type_list)):
            if cell_type_list[j] == unique_types[i]:
                cell_type_markers[unique_types[i]].append(marker_list[j])
                j += 1
                
    return(cell_type_markers)


### load Cusanovich peak/promoter intersections
### need to also set up function to create the file from scratch instead of taking it from Cusanovich 
### (I'm pretty sure I anyway already have a script for that)

def identify_cluster(adata, cell_type, cell_type_markers, peak_promoter_file, gene_name_pos=5, path='', n_peaks_per_cluster=1000):
    """
    Use markers of a given cell type to plot peak openness for peaks in promoters of the given markers
    Input cell type, cell type markers, peak promoter intersections 
    
    Paramters
    ---------
    cell_type: str, cell type that is to be investigated (must be the same as in the dictionary)
    cell_type_markers: dict, output of load_markers
    peak_promoter_file: tab-separated file that stores information about peak/promoter intersections and gene names for these promoters
    gene_name_pose: int, indicates the column in the peak_promoter_file where the gene name is stored
    path: str, path to the peak_promoter_file
    n_peaks_per_cluster: int, number of peaks per louvain cluster that should be searched for matches with the markers
    
    Return
    ------
    umap depicting peak openness for promoters of marker genes for a given cell type
    """
    if "omic" in adata.uns.keys():
        if adata.uns['omics'] == 'ATAC':
            return(identify_cluster_ATAC(adata, cell_type, cell_type_markers,
                peak_promoter_file, gene_name_pos, path, n_peaks_per_cluster))
        else:
            message = 'the current implementation of this function only allow ATAC-seq data as input'
            warnings.warn(message)
    else:
        message = """the current implementation of this function only allow ATAC-seq data as input.\n
        Additionally, annData uns['omics'] is not specified."""
        warnings.warn(message) 

    return()



def identify_cluster_ATAC(adata, cell_type, cell_type_markers, peak_promoter_file, gene_name_pos=5, path='', n_peaks_per_cluster=1000):
    
    peak_name = []
    gene_name = []
    with open(path+peak_promoter_file) as f:
        head = f.readline()
        for line in f:
            line = line.split("\t")
            peak_name.append(line[0])
            gene_name.append(line[gene_name_pos].upper())
            
    matching_peaks = []
    matching_gene_name = []
    for i in range(0,len(peak_name)):
        # check if peak is in the adata object and if the according gene name is a marker
        if peak_name[i] in list(adata.var_names) and gene_name[i] in cell_type_markers[cell_type]:
            matching_peaks.append(peak_name[i])
            matching_gene_name.append(gene_name[i])
    # get the cell type names per matching peak
    cell_type_peaks = {}
    cell_type_peaks[cell_type] = []
    for j in range(0,len(matching_peaks)):
        if matching_gene_name[j] in cell_type_markers[cell_type]:
            cell_type_peaks[cell_type].append(matching_peaks[j]+"_"+matching_gene_name[j])
            
    # make the peaks unique
    cell_type_peaks[cell_type] = list(set(cell_type_peaks[cell_type]))
        
        
    cell_type_peak = []
    cell_type_gene = []
    for elem in cell_type_peaks[cell_type]:
        ctype = elem.split("_")
        cell_type_peak.append(ctype[0]+"_"+ctype[1]+"_"+ctype[2])
        cell_type_gene.append(ctype[3])
    
    sc.tl.rank_genes_groups(adata, groupby="louvain", n_genes=n_peaks_per_cluster)

    # top 100 to 1000 ranked peaks per louvain group
    ATAC_ranking_dict = {}
    for group in range(0, len(set(list(adata.obs["louvain"])))):
        ATAC_ranking_dict[str(group)] = []
        for idx in range(0,n_peaks_per_cluster):
            ATAC_ranking_dict[str(group)].append(list(adata.uns["rank_genes_groups"]['names'])[idx][group])
        
    sig_markers = []
    sig_peak = []
    for group in range(0, len(set(list(adata.obs["louvain"])))):
        for elem in ATAC_ranking_dict[str(group)]:
            if elem in cell_type_peak:
                sig_markers.append(elem+":"+str(group))
                sig_peak.append(elem)                
            
    if len(sig_markers) > 200:
        print("WARNING: found more than 200 marker peaks for your chosen cell type. Consider choosing a smaller number for n_peaks_per_cluster to reduce running time and get more significant results.")
    
    group_contribution = []
    for element in sig_markers:
        group = element.split(":")[1]
        group_contribution.append(group)

    most_prominent = max(set(group_contribution), key = group_contribution.count) 
    
    print("Differentially open peaks in promoters for known",cell_type,"marker genes are most prominently found in louvain group",most_prominent)
    
    # only take unique peaks (peaks can be among top markers in more than one group)
    unique_sig_peaks = []
    sig_gene_names = []
    for idx,name in enumerate(cell_type_peak):
        if name in sig_peak:
            sig_gene_names.append(cell_type_gene[idx])
            unique_sig_peaks.append(cell_type_peak[idx])
        
    return(sc.pl.umap(adata, color=unique_sig_peaks, title=sig_gene_names))   