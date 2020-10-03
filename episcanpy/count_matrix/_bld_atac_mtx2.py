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
        
def parallel_counting(idx_parts, window_list, barcodes, index):
    tbx = pysam.TabixFile(bed_file)
    mtx = lil_matrix((len(barcodes), len(window_list)), dtype=np.uint16)
    print("In Index: ", index, " PID: ", os.getpid(), ", loading tbx file")
    print(idx_parts[index])
    for i in idx_parts[index]:
        for row in tbx.fetch(window_list[i][0], window_list[i][1], window_list[i][2], parser=pysam.asTuple()):
            mtx[barcodes.index(str(row).split('\t')[3].split(':')[0]), i] += 1
    allmtx[0] = allmtx[0] + mtx
    ss=mtx.sum().sum()
    print("In Index: ", index, " PID: ", os.getpid(),", all value = ",ss)
                
# the codes are from EpiScanpy but adapt for these input files by using pybedtools
def bld_mtx_fly(bed_file, annotation, chrom, csv_file=None, genome=None, thread=1, save=False):
    """
    Building count matrix on the fly.
    Expected running time for 10k cells X 100k features on a personal computer ~65min
    Does not count pcr duplicate.
    A tbi file with the same filename as the provided bed_file must be present in the same directory as the tsv file
    Parameters
    ----------
    bed_file : name of the file containing the multiplexed reads (.tsv or .tsv.gz)
    annotation : loaded set of features to create the feature matrix from
    chrom : chromosome
    csv_file : default is None -
    genome : default is None - specify if you want to extract a specific genome assembly
    save : default is False - supply a file path as str to save generated AnnData object
    Output
    ------
    AnnData object (also saved as h5ad if save argument is specified)
    """
    start = time.time()
    manager = Manager()
    global allmtx
    global bed_file
    allmtx = manager.dict()
    intervaltime = time.time()
    print("Time point, loading barcodess " + str(intervaltime-start) + " sec")
    
    key = 'chr' + chrom
    print("Chr: ", key)
    #barcodes = pd.DataFrame(columns=[0,1,2,3,4])
    barcodes = set()
    #idx = 0
    with gzip.open(bed_file,'rt') as f:
        for line in f:
            ar = line.strip().split('\t')
            if ar[0] == key:
                barcodes.add(ar[3])
                #idx += 1
    #print(barcodes)
    barcodes = sorted(list(barcodes))
    #print(barcodes)
    #intervaltime = time.time()
    #print("Time point, sorting barcodes " + str(intervaltime-start) + " sec")
    #barcodes = sorted(list(set([x.split(':')[0] for x in barcodes.loc[:, 3].unique().tolist()])))
    #print(barcodes)
    #barcodes = sorted(list(set([x.split(':')[0] for x in pd.read_csv(bed_file, sep='\t', header=None).loc[:, 3].unique().tolist()])))
    
    intervaltime = time.time()
    print("Time point, done sorting barcodes " + str(intervaltime-start) + " sec")
    print("barcodes:\t", len(barcodes))
    
    # barcodes
    #nb_barcodes = len(barcodes)
    #dict_barcodes = {barcodes[i]: i for i in range(0, len(barcodes))}
    
    # Load tabix
    intervaltime = time.time()
    print("Time point, loading index file " + str(intervaltime-start) + " sec")
    #tbx = pysam.TabixFile(bed_file)
    intervaltime = time.time()
    print("Time point, done loading index file and creating window_list " + str(intervaltime-start) + " sec")
    
    # format annotations
    #chromosome_list =['chr'+x for x in annotation.keys()]  # Make a list of chromosomes
    window_list = []
    #print(annotation.keys())
    if genome:
        window_list += [["".join([genome, '_chr', str(chrom)]), int(n[0]), int(n[1])] for n in annotation[chrom]]
    else:
        window_list += [["".join(['chr', str(chrom)]), int(n[0]), int(n[1])] for n in annotation[chrom]]

    intervaltime = time.time()
    print("Time point, after creating window_list " + str(intervaltime-start) + " sec")
    
    if genome:
        for chrom in sorted(annotation.keys()):
            window_list += [["".join([genome, '_chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]
    else:
        for chrom in sorted(annotation.keys()):
            window_list += [["".join(['chr', chrom]), int(n[0]), int(n[1])] for n in annotation[chrom]]

    print('building count matrix')
    #mtx = lil_matrix((len(barcodes), len(window_list)), dtype=np.uint16)
    #print(mtx.shape)
    
    intervaltime = time.time()
    print("Time point, after creating empty matrix " + str(intervaltime-start) + " sec")
    
    
    #Parallel(n_jobs=num_cores, prefer="threads")(delayed (parallel_counting)(i) for i in range(len(window_list)))
    idx_parts = chunkIt(range(len(window_list)),thread)
    print("All parts for threading")
    print(idx_parts)
    allmtx[0] = lil_matrix((len(barcodes), len(window_list)), dtype=np.uint16)
    p = Pool(thread)
    func = partial(parallel_counting, idx_parts, window_list, barcodes)
    p.map(func, range(len(idx_parts)))
    p.close()
    p.join()
    
    #mtx = allmtx[0]
    #print("All mtx ",len(allmtx))
    #mtx = lil_matrix((len(barcodes), len(window_list)), dtype=np.uint16)
    #for i in allmtx:
    #    mtx = mtx + allmtx[i]
    ss=allmtx[0].sum().sum()
    print("All value = ",ss)
    #with Pool(THREAD) as p: # the argument in Pool is the number of processes
    #    p.map(parallel_counting(idx_parts), range(len(idx_parts)))
    
    # for i in range(len(window_list)):
    #     for row in tbx.fetch(window_list[i][0], window_list[i][1], window_list[i][2], parser=pysam.asTuple()):
    #         mtx[barcodes.index(str(row).split('\t')[3].split(':')[0]), i] += 1

#    for i, tmp_feat in enumerate(tqdm(window_list)):
#        intervaltime = time.time()
#        print("Time point, in window_list loop of ", i," and  ", tmp_feat, "   " + str(intervaltime-start) + " sec")
#        for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
#            mtx[barcodes.index(str(row).split('\t')[3].split(':')[0]), i] += 1
    
    intervaltime = time.time()
    print("Time point, after creating a count matrix " + str(intervaltime-start) + " sec")
    
    print('building AnnData object')
    allmtx[0] = ad.AnnData(allmtx[0].tocsr(),
                     obs=pd.DataFrame(index=barcodes),
                     var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in window_list]))

    intervaltime = time.time()
    print("Time point, after saving H5AD " + str(intervaltime-start) + " sec")
    
    if csv_file:
        print('filtering barcodes')
        df = pd.read_csv(csv_file)
        if genome == 'mm10':
            df_filtered = df[(df.is_mm10_cell_barcode == 1)]
        elif genome == 'hg19':
            df_filtered = df[(df.is_hg19_cell_barcode == 1)]
        else:
            df_filtered = df[(df.is_cell_barcode == 1)]

        barcodes = set(df_filtered.barcode.tolist())
        allmtx[0] = allmtx[0][[i in barcodes for i in allmtx[0].obs.index]].copy()

    if save:
        allmtx[0].write(save)

    return(allmtx[0])