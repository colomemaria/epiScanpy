import pysam
import numpy as np
import anndata as ad
import pandas as pd

MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']

def bld_atac_mtx(list_bam_files, loaded_feat, output_file='', path='', writing_option='a', header=None, chromosomes=MOUSE):
    """
    Differently to the methylation function to build count matrix, you can only write one
    set of data at any given time. I need to change this to allow to build multiple matrices at the same time. 
    
    chromose is either MOUSE, HUMAN or a list of chromosomes. 
    MOUSE=['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
    HUMAN=['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']
    
    """
        
    if output_file=='':
        output_file='std_output_ct_mtx.txt'
    # open file to write
    output_file = open(path+output_file, writing_option)
    
    # write header if specified
    if header != None:
        output_file.write('sample_name\t')
        for feature in header:
            output_file.write(feature)
            output_file.write('\t')
        output_file.write('\n')
    
    # start going through the bam files
    for name_file in list_bam_files[0:]:
        ## important variables for output
        index_feat = {key: 0 for key in chromosomes}    
        val_feat = {key: [0 for x in range(len(loaded_feat[key]))] for key in chromosomes}
    
        ## PART 1 read the bam file
        keep_lines = []
        samfile = pysam.AlignmentFile(path+name_file, "rb")
        for read in samfile.fetch(until_eof=True):
            line = str(read).split('\t')
            if line[2] in chromosomes:
                keep_lines.append(line[2:4])
            ### print -- output
        print(name_file, len(keep_lines), 'mapped reads')
        samfile.close()
        
        ## PART2 reads that fall into 
        for element in keep_lines:
            ## 2 things per line:
            chrom = element[0]
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
        # write down the result of the cell
        output_file.write(name_file)
        output_file.write('\t')
        for chrom in chromosomes:
            output_file.write('\t'.join([str(p) for p in val_feat[chrom]]))
            output_file.write('\t')
        output_file.write('\n')

    output_file.close()
    
#convert the matrix into a sparse matrix:

def save_sparse_mtx(initial_matrix, output_file='.h5ad', path='', omic='ATAC'):
    """
    I need to add the methylation conversion to sparse matrix. Right now there is no need for
    the omic argument because I haven't added the conversion for methylation. It only deal with "dropout omics"
    like ChIPseq, ATACseq or RNAseq
    """
    head = None
    data = []
    cell_names = []
    # reading the non sparse file
    with open(path+initial_matrix) as f:
        first_line = f.readline()
        first_line = first_line[:-3].split('\t')
        if first_line[0] == 'sample_name':
            head = first_line
        else:
            cell_names.append(first_line[0])
            data = [[int(l) for l in first_line[1:]]]
        file = f.readlines()
    
    for line in file:
        line = line[:-3].split('\t')
        cell_names.append(line[0])
        data.append([int(l) for l in line[1:]])

    # convert into an AnnData object
    if head != None:
        adata = ad.AnnData(np.array(data), obs=pd.DataFrame(index=cell_names), var=pd.DataFrame(index=head[1:]))
    else:
        adata = ad.AnnData(np.array(data), obs=pd.DataFrame(index=cell_names))
        
    # writing the file as h5ad --> sparse matrix with minimum annotations
    if output_file=='.h5ad':
        output_file = "".join([initial_matrix.split('.')[0], output_file])
        
    adata.write(output_file)
    
