#import pysam
import bamnostic as bs
import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse  import  csc_matrix

MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']

def bld_atac_mtx(list_bam_files, loaded_feat, output_file_name=None,
    path=None, writing_option='a', header=None, mode='rb',
    check_sq=True, chromosomes=HUMAN):
    """
    Build a count matrix one set of features at a time. It is specific of ATAC-seq data.
    It curently do not write down a sparse matrix. It writes down a regular count matrix
    as a text file. 
    
    Parameters
    ----------

    list_bam_files: input must be a list of bam file names. One for each cell to 
        build the count matrix for

    loaded_feat: the features for which you want to build the count matrix
        
    output_file_name: name of the output file. The count matrix that will be written
        down in the current directory. If this parameter is not specified, 
        the output count amtrix will be named 'std_output_ct_mtx.txt'

    path: path where to find the input file. The output file will be written down
    in your current directory, it is not affected by this parameter.

    writing_option: standard writing options for the output file. 'a' or 'w'
        'a' to append to an already existing matrix. 'w' to overwrite any 
        previously exisiting matrix. 
        default: 'a'

    header: if you want to write down the feature name specify this argument.
        Input must be a list.

    mode: bamnostic argument 'r' or 'w' for read and write 'b' and 's' for bam or sam
        if only 'r' is specified, bamnostic will try to determine if the input is 
        either a bam or sam file.

    check_sq: bamnostic argument. when reading, check if SQ entries are present in header

    chromosomes: chromosomes of the species you are considering. default value
        is the human genome (not including mitochondrial genome).
        HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
                '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
        MOUSE = '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']

    Return
    ------
    It does not return any object. The function write down the desired count
    matrix in a txt file

    """
        
    if output_file_name==None:
        output_file_name='std_output_ct_mtx.txt'


    if path==None:
        path=''
    
    # open file to write
    output_file = open(path+output_file_name, writing_option)
    # write header if specified
    if header != None:
        output_file.write('sample_name\t')
        for feature in header:
            output_file.write(feature)
            output_file.write('\t')
        output_file.write('\n')
    # close file to write   
    output_file.close()

    # start going through the bam files
    for name_file in list_bam_files[0:]:

        ## important variables for output
        index_feat = {key: 0 for key in chromosomes}    
        val_feat = {key: [0 for x in range(len(loaded_feat[key]))] for key in chromosomes}
    
        ## PART 1 read the bam file
        keep_lines = []
        #samfile = bs.AlignmentFile(path+output_file_name, mode="rb", check_sq=False)
        samfile = bs.AlignmentFile(path+name_file, mode="rb", check_sq=False)
        #for read in samfile.fetch(until_eof=True):
        for read in samfile:
            line = str(read).split('\t')
            if line[2][3:] in chromosomes:
                keep_lines.append(line[2:4])
            ### print -- output
        print(name_file, len(keep_lines), 'mapped reads')
        samfile.close()
        
        ## PART2 reads that fall into 
        for element in keep_lines:
            ## 2 things per line:
            chrom = element[0][3:]
            read_pos = int(element[1])
            max_value_index = len(loaded_feat[chrom])
            ## I want to check if the read map to a feature in the same chrom
            pointer_feat_pos = index_feat[chrom]
            for feat_pos in loaded_feat[chrom][pointer_feat_pos:]:
                pointer_feat_pos += 1
                # Go through all features that are smaller than the read position
                if read_pos > feat_pos[1]:
                    continue
                # if read_pos fall in a feature
                elif read_pos > feat_pos[0]:
                    # Update the pointer for the next read if the pointer isn't out of range
                    if pointer_feat_pos < max_value_index:
                        index_feat[chrom] = pointer_feat_pos
                        val_feat[chrom][pointer_feat_pos] += 1
                    else:
                        index_feat[chrom] = max_value_index
                    # Check the following features without updating the pointer. 
                    break
                else:
                    break
     
            for feat_pos in loaded_feat[chrom][pointer_feat_pos:]:
                # +1 if reads fall into more than one feature
                if feat_pos[0] < read_pos:
                    val_feat[chrom][pointer_feat_pos] += 1
                    pointer_feat_pos += 1
                # if read_pos > start position of the new feature break
                elif read_pos < feat_pos[0]:
                    break
                else:
                    print('error')
                    break
                    
        # PART 3
        # open
        output_file = open(path+output_file_name, 'a')
        # write down the result of the cell
        output_file.write(name_file)
        output_file.write('\t')
        for chrom in chromosomes:
            output_file.write('\t'.join([str(p) for p in val_feat[chrom]]))
            output_file.write('\t')
        output_file.write('\n')
        #close
        output_file.close()
    
def read_mtx_bed(file_name, path='', omic='ATAC'):
    """
    read this specific matrix format. It is the standard output of bedtools when you merge bam files.
    """
    peak_name = []
    cell_matrix = [] 
    with open(path+file_name) as f:
        head = f.readline().split('\t')
        head[len(head)-1] = head[len(head)-1].split("\n")[0]
        for line in f:
            line = line.split('\t')
            line[len(line)-1] = line[len(line)-1].split("\n")[0]
            peak_name.append(line[3]) # for some reason it has rownames
            cell_matrix.append([int(x) for x in line[4:]])
    cell_names = head[4:]
    cell_matrix=np.matrix(cell_matrix)
    cell_matrix = cell_matrix.transpose()
    adata = ad.AnnData(cell_matrix,
                   obs=pd.DataFrame(index=cell_names),
                   var=pd.DataFrame(index=peak_name))
    if omic != None:
        adata.uns['omic'] = omic
    return(adata)


def save_sparse_mtx(initial_matrix, output_file='.h5ad', path='', omic='ATAC', bed=False, save=True):
    """
    Convert regular atac matrix into a sparse Anndata:

    Parameters
    ----------

    initial_matrix: initial dense count matrix to load and convert into a sparse matrix 

    output_file: name of the output file for the AnnData object.
    Default output is the name of the input file with .h5ad extension

    path: path to the input count matrix. The AnnData object is written in the current directory, 
        not the location specified in path.

    omic: 'ATAC', 'RNA' or 'methylation' are the 3 currently recognised omics in epiScanpy. 
        However, other omic name can be accepted but are not yet recognised in other functions.
        default: 'ATAC'

    bed: boolean. If True it consider another input format (bedtools output format for count matrices)

    save: boolean. If True, the sparse matrix is saved as h5ad file. Otherwise it is simply return.

    Return
    ------

    It returns the loaded matrix as an AnnData object.
    """
    head = None
    data = []
    cell_names = []

    # choice between 2 different input count matrix formats
    if bed == True:
        adata = read_mtx_bed(file_name, path, omic)
    else:
        # reading the non sparse file
        with open(path+initial_matrix) as f:
            first_line = f.readline()
            first_line = first_line[:-3].split('\t')
            if first_line[0] == 'sample_name':
                head = first_line[:-1]
            else:
                cell_names.append(first_line[0])
                data = [[int(l) for l in first_line[1:-1]]]
            file = f.readlines()
    
        for line in file:
            line = line[:-3].split('\t')
            cell_names.append(line[0])
            data.append([int(l) for l in line[1:-1]])
            

        # convert into an AnnData object
        if head != None:
            adata = ad.AnnData(csc_matrix(data), obs=pd.DataFrame(index=cell_names), var=pd.DataFrame(index=head[1:]))
        else:
            adata = ad.AnnData(csc_matrix(data), obs=pd.DataFrame(index=cell_names))
        
        if omic != None:
            adata.uns['omic'] = omic

    # writing the file as h5ad --> sparse matrix with minimum annotations
    if save ==True:
        if output_file=='.h5ad':
            output_file = "".join([initial_matrix.split('.')[0], output_file])
        
        adata.write(output_file)

    return(adata)