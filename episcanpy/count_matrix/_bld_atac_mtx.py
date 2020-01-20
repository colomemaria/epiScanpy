import numpy as np
import anndata as ad
import pandas as pd
import gzip
import pysam
from scipy.sparse import csc_matrix

def bld_mtx_fly(tsv_file, annotation, csv_file=None, genome=None, DATADIR='', save=False):
    """
    Building count matrix on the fly. 
    Expected running time for 10k cells X 100k features on a personal computer ~65min
    Does not count pcr duplicate. 
    Input must be tsv file with tbi
    
    Parameters
    ----------
    
    tsv_file : name of the file containing the multiplexed reads.
    
    annotation : loaded set of features to 
    
    csv_file : default is None - 
    
    genome : default is None - specify if you want to extract a specific genome assembly
        (if )
    
    DATADIR : default is '' - 
    
    save : default is False - 
    
    Output
    ------
    
    Adata object  (saved as h5ad)
    
    """
    
    if csv_file != None:
        print('loading and filtering the barcodes')
        df = pd.read_csv(DATADIR+csv_file)
        if genome=='mm10':
            # Filter by is__cell_barcode == 1
            df_filtered = df[(df.is_mm10_cell_barcode == 1)]
        elif genome=='hg19':
            df_filtered = df[(df.is_hg19_cell_barcode == 1)]
        else:
            df_filtered = df[(df.is__cell_barcode == 1)]
        
        df_sorted = df_filtered.sort_values('passed_filters', inplace=False, ascending=False) 
        # Get the barcodes list as an array
        barcodes = df_sorted.barcode.values
        
    else:
        print('loading the barcodes')
        with gzip.open(tsv_file, 'rb') as f_in:
            with open(tsv_file.rstrip('.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        df = pd.read_csv(tsv_file.rstrip('.gz'), sep='\t', header=None)
        barcodes = list(sorted(set(df.loc[:,3].tolist())))
        
    # barcodes
    nb_barcodes = len(barcodes)
    dict_barcodes = {barcodes[i] : i for i in range(0, len(barcodes))}
    
    # Load tabix
    tbx = pysam.TabixFile(DATADIR+tsv_file)
    
    print('building the count matrix')
    mtx = []

    for tmp_feat in annotation: 
        vector = [0]*nb_barcodes
        for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
            line = str(row).split('\t')[-2]
            if line in barcodes:
                vector[dict_barcodes[line]] += 1 
        
        mtx.append(vector)
    print('making the AnnData')
    
    mtx = ad.AnnData(csc_matrix(np.array(mtx)).transpose(),
                     obs=pd.DataFrame(index=barcodes),
                     var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in annotation]))
    if save!=False:
        mtx.write(save)
    
    return(mtx)

def bld_mtx_fly(tsv_file, tbi_file, annotation, csv_file=None, genome=None, DATADIR='', save=False):
    """
    Building count matrix on the fly. 
    Expected running time for 10k cells X 100k features on a personal computer ~65min
    Does not count pcr duplicate. 
    Input must be tsv file with tbi
    
    Parameters
    ----------
    
    tsv_file : name of the file containing the multiplexed reads.
    
    annotation : loaded set of features to 
    
    csv_file : default is None - 
    
    genome : default is None - specify if you want to extract a specific genome assembly
        (if )
    
    DATADIR : default is '' - 
    
    save : default is False - 
    
    Output
    ------
    
    Adata object  (saved as h5ad)
    
    """
    if csv_file != None:
        print('loading and filtering the barcodes')
        df = pd.read_csv(DATADIR+csv_file)
        if genome=='mm10':
            # Filter by is__cell_barcode == 1
            df_filtered = df[(df.is_mm10_cell_barcode == 1)]
        elif genome=='hg19':
            df_filtered = df[(df.is_hg19_cell_barcode == 1)]
        else:
            df_filtered = df[(df.is__cell_barcode == 1)]
        
        df_sorted = df_filtered.sort_values('passed_filters', inplace=False, ascending=False) 
        # Get the barcodes list as an array
        barcodes = df_sorted.barcode.values
        
    else:
        print('loading the barcodes')
        with gzip.open(tsv_file, 'rb') as f_in:
            with open(tsv_file.rstrip('.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        df = pd.read_csv(tsv_file.rstrip('.gz'), sep='\t', header=None)
        barcodes = list(sorted(set(df.loc[:,3].tolist())))
        
    # barcodes
    nb_barcodes = len(barcodes)
    dict_barcodes = {barcodes[i] : i for i in range(0, len(barcodes))}
    
    # Load tabix
    tbx = pysam.TabixFile(DATADIR+tsv_file)
    
    # format annotations
    window_list = []
    for chrom in sorted(annotation.keys()):
        window_list += [[chrom, int(n[0]), int(n[1])] for n in annotation[chrom]]
    
    print('building the count matrix')
    mtx = []
    for tmp_feat in window_list: 
        vector = [0]*nb_barcodes
        for row in tbx.fetch(tmp_feat[0], tmp_feat[1], tmp_feat[2], parser=pysam.asTuple()):
        #for row in tbx.fetch(start=tmp_feat[0], end=tmp_feat[1], region=tmp_feat[2], parser=pysam.asTuple()):
            line = str(row).split('\t')[-2]
            if line in barcodes:
                vector[dict_barcodes[line]] += 1 
        
        mtx.append(vector)
    print('making the AnnData')
    mtx = csc_matrix(np.array(mtx))
    mtx = ad.AnnData(mtx.transpose(),
                     obs=pd.DataFrame(index=barcodes),
                     var=pd.DataFrame(index=['_'.join([str(p) for p in n]) for n in window_list]))
    if save!=False:
        mtx.write(save)
    
    return(mtx)