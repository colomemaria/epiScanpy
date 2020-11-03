import time
import pysam
import argparse
import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import lil_matrix
from multiprocessing import Pool, Manager, Value
import gzip
import os
from functools import partial

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out
        
def parallel_counting(bed_file,allmtx, idx_parts, feature_list, cell_ids, cell_id_col, readname_sep, index):
    tbx = pysam.TabixFile(bed_file)
    mtx = lil_matrix((len(cell_ids), len(feature_list)), dtype=np.uint16)
    #print("In Index: ", index, " PID: ", os.getpid(), ", loading tbx file")
    #print(idx_parts[index])
    for i in idx_parts[index]:
        for row in tbx.fetch(feature_list[i][0], feature_list[i][1], feature_list[i][2], parser=pysam.asTuple()):
            # print(str(row))
            if readname_sep:
                mtx[cell_ids.index(str(row).split()[3].split(readname_sep)[cell_id_col]), i] += 1
            else:
                mtx[cell_ids.index(str(row).split()[3]), i] += 1
    allmtx[0] = allmtx[0] + mtx
    ss=mtx.sum().sum()
    tbx.close()
    #print("In Index: ", index, " PID: ", os.getpid(),", all value = ",ss)
                

def bld_mtx_bed_per_chr(bed_file, feature_region, chrom, cell_id_col = 0, readname_sep = ':', thread=1, save=False):
    """
    Building count matrix on the fly.
    Expected running time for 10k cells X 100k features on a personal computer ~65min
    Does not count pcr duplicate.
    A tbi file with the same filename as the provided bed_file must be present in the same directory as the tsv file
    Parameters
    ----------
    bed_file : a path to the file containing the multiplexed reads (.tsv or .tsv.gz)
    feature_region: a dictionary containing the feature regions built by the function make_windows or a path to the file containing the predefined regions for all feattures. The file must contain 3 columns; chromosome, start position, and end position, separated by tab or space. chromosome should be in the format, as "chrom1", "chrom2", ..., "chromY", "chromM".
    chrom : chromosome
    cell_id_col : a column index of cell identity on the read name. For example, 'AAATAA:NT500:65:EXP01', you know that the cell identity is at the first column, and it is separated by ':', you can set cell_id_col = 0. If 'AAATAA:NT500:65:EXP01' is a cell identiy, you can set cell_id_col = 0 and readname_sep = None.
    readname_sep : : a separator on the read name. For example, 'AAATAA:NT500:65:EXP01', you know that the cell identity is at the first column, and it is separated by ':', you can set readname_sep = ':'. However, if 'AAATAA:NT500:65:EXP01' is a cell identiy, you can set readname_sep = None.
    thread : a number of parallel threads to run.
    save : default is False - supply a file path as str to save generated AnnData object
    Output
    ------
    AnnData object (also saved as h5ad if save argument is specified)
    """
    start = time.time()
    manager = Manager()
    global allmtx
    allmtx = manager.dict()
    
    key = 'chr' + chrom
    # print("chr: ", key)
    
    feature = {}
    if not(isinstance(feature_region, dict)): 
        feature[chrom] = []
        file = open(feature_region)
        for line in file:
            ar = line.strip().split()
            if ar[0] == ('chr' + chrom):
                feature[chrom].append([int(ar[1]), int(ar[2])])
        file.close()
    else:
        feature = feature_region
        del(feature_region)
    
    read_names = set()
    with gzip.open(bed_file,'rt') as f:
        for line in f:
            ar = line.strip().split()
            if ar[0] == key:
                read_names.add(ar[3])
    
    if readname_sep:
        cell_ids = sorted(list(set([n.split(readname_sep)[cell_id_col] for n in read_names])))
    else:
        cell_ids = sorted(list(set(read_names)))
    
    # intervaltime = time.time()
    # print("Time point, after fetching cell ID " + str(intervaltime-start) + " sec")
    # print("Number of cells:\t", len(cell_ids))
    
    feature_list = []
    feature_list += [["".join(['chr', str(chrom)]), int(n[0]), int(n[1])] for n in feature[chrom]]
    
    idx_parts = chunkIt(range(len(feature_list)),thread)
    # print("All parts for threading")
    # print(idx_parts)
    allmtx[0] = lil_matrix((len(cell_ids), len(feature_list)), dtype=np.uint16)
    p = Pool(thread)
    func = partial(parallel_counting, bed_file, allmtx, idx_parts, feature_list, cell_ids, cell_id_col, readname_sep)
    p.map(func, range(len(idx_parts)))
    p.close()
    p.join()
    
    intervaltime = time.time()
    # print("Time point, after creating a count matrix " + str(intervaltime-start) + " sec")
    ss=allmtx[0].sum().sum()
    # print("All valid reads = ",ss)
    allmtx[0] = ad.AnnData(allmtx[0].tocsr(),
                     obs=pd.DataFrame(index=cell_ids),
                     var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in feature_list]))
    
    if save:
        allmtx[0].write(save)
        
    intervaltime = time.time()
    print("Chromosome " + chrom + ": " + str(intervaltime-start) + " sec")
    
    return(allmtx[0])

def bld_mtx_bed(bed_file, feature_region = None, cell_id_col = 0, readname_sep = ':', chromosomes = 'human', thread=1, save = False, bin_size = 5000):
    """
    Building count matrix from BED file.
    Does not count pcr duplicate.
    A tbi file with the same filename as the provided bed_file must be present in the same directory as the BED file
    Parameters
    ----------
    bed_file : a path to the file containing the reads (.bed or .bed.gz)
    feature_region: a dictionary containing the feature regions built by the function make_windows or a path to the file containing the predefined regions for all feattures. Default value is None. 
    If chromosomes is set as 'human' or 'mouse', it will be unused. The file must contain 3 columns; chromosome, start position, and end position, separated by tab or space. chromosome should be in the format, as "chrom1", "chrom2", ..., "chromY", "chromM".
    cell_id_col : a column index of cell identity on the read name. For example, 'AAATAA:NT500:65:EXP01', you know that the cell identity is at the first column, and it is separated by ':', you can set cell_id_col = 0. If 'AAATAA:NT500:65:EXP01' is a cell identiy, you can set cell_id_col = 0 and readname_sep = None.
    readname_sep : : a separator on the read name. For example, 'AAATAA:NT500:65:EXP01', you know that the cell identity is at the first column, and it is separated by ':', you can set readname_sep = ':'. However, if 'AAATAA:NT500:65:EXP01' is a cell identiy, you can set readname_sep = None.
    chromosomes : chromosomes of the species you are considering. It can be set as 'human' and 'mouse' and the default value is 'human'.
    'human' refers to ['1', '2', '3', ... , '22', 'X', 'Y']
    'mouse' refers to ['1', '2', '3', ... ', '19', 'X', 'Y']
    For the other species, chromosomes must be defined as, for example,
    chromosomes = ['1', '2', '3', ... ,'X', 'Y']
    For mitochondrial, it should be defined as 'M' or 'MT' and it must be matched with the column chromosome in feature_region.
    thread : a number of parallel threads to run.
    save : default is False - supply a file path as str to save generated AnnData object
    bin_size : the size of non-overlapped windows. Default value is 5000. It is used when feature_region is None, or chromosomes is set as 'human' or 'mouse'.
    Output
    ------
    AnnData object (also saved as h5ad if save argument is specified)
    """
    start = time.time()
    MOUSE = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
            '11', '12', '13', '14', '15', '16', '17', '18', '19','X', 'Y']
    HUMAN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','X', 'Y']
    if not(isinstance(chromosomes, list)):
        if chromosomes == 'human':
            chromosomes = HUMAN
            feature_region = make_windows(size=bin_size, chromosomes='human')
        if chromosomes == 'mouse':
            chromosomes = MOUSE
            feature_region = make_windows(size=bin_size, chromosomes='mouse')
    if not(feature_region):
        print("feature_region must not be None when chromosomes is not 'human' or 'mouse'. Please set feature_region to a path of feature file or generate from the funtion make_windows.")
    adata_outter = None
    for ch in chromosomes:
        tmp = bld_mtx_bed_per_chr(bed_file = bed_file, 
                                feature_region = feature_region,
                                chrom = ch,
                                cell_id_col = cell_id_col,
                                readname_sep = readname_sep,
                                thread = thread,
                                save = False)
        print("Chromosome " + ch)
        # print(tmp)
        if len(tmp.var) != 0:
            if adata_outter is None:
                adata_outter = tmp.copy()
            else:
                adata_outter = ad.concat([adata_outter, tmp], axis=1, join = "outer")
        else:
            print("Chromosome " + ch + " contains no data")
    print("All data contains:")
    print(adata_outter)
    if save:
        adata_outter.write(save)
        
    intervaltime = time.time()
    print("Total time is " + str(intervaltime-start) + " sec")
    return(adata_outter)

