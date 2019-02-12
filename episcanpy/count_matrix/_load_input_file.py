import numpy as np
from numpy import NaN
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True)


# chromosomes for 2 principal species. If you work with another genome
# the chromosomes will have to be specified
HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2',
         '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2',
         '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
# my favorite dataset ;-)
HEAD_Ecker = ['chr', 'pos', 'strand', 'mc_class', 'mc_count', 'total', 'methylated\n']


def read_cyt_summary(sample_name, meth_type, head, path, chromosome):
    """
    Read file from which you want to extract the methylation level and
    (assuming it is like the Ecker/Methylpy format) extract the number of
    methylated read and the total number of read for the cytosines covered and
    in the right genomic context (CG or CH)
    Parameters
    ----------
    sample_name:
        name of the file to read to extract key information.
    meth_type:
        CG, CH or not specified
    head: if there is header that you don't want to read. An annotation in the
        file you are reading. The default value is the Methylpy/Ecker header
    path: path of the to access the file of the sample you want to read. 
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes)
    """
    # open the sample file
    with open(path+sample_name) as sample:
        
        # initialize the variables
        reduced_cyt = {} # to store cyt for each chrom (intermed output)
        for c in chromosome:
            reduced_cyt[c] = []

        # specify the cytosine context you want to consider:
        if meth_type == 'CG':
            extract_CG(sample, reduced_cyt, head)
        elif meth_type == 'CH':
            extract_CH(sample, reduced_cyt, head)
        elif meth_type == '':
            for line in sample:
                line = line.split("\t")
                if line != head:
                    reduced_cyt[line[0]].append([int(line[1]),int(line[-3]), int(line[-2])])
        else:
            print ("error in the choice of methylation context")
            return()
    return(reduced_cyt)
