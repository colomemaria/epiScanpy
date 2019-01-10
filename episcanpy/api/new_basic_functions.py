#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 13:16:53 2018
@author: Anna Danese

We gather here the basic fucntions necessaries to process files to get features, 
and extract methylation level from a file. 
And produce some summary statistics. 

PART 1: prepare features and windows to measure methylation levels 
PART 2: process file
PART 3: summary statistics
PART 4: additional processing of data (methylation levels and matrices)

"""

################################ PART 0 #######################################

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

###############################################################################
################################ PART 1 #######################################
###############################################################################

# it assume that you take a bed file as impute
# it save the features position and name like [start, end, name]
# in a dictionary that has chromosomes as keys. 
def load_features(file_features, chromosomes=HUMAN, path="", sort=False):
    """
    The function load features is here to transform a bed file into a usable 
    set of units to measure methylation levels. 
    It has to be a bed-like file. You also need to specify the chromosomes you 
    use as a list of characteres like ['1', '7', '8', '9', 'X', 'Y', 'M'].
    The chromosome list you give as input can be not ordered.
    If you don't specify the chromosomes, the default is the human genome
    (including X, Y and mitochondrial DNA).
    THE BED FILE need to be sorted !! (Maybe I should add an option to sort the
    file).
    
    The output is a dictionary where the keys are chromosomes and the value is
    a list containing [start, end, name] for every feature extracted.
    
    Parameters
    ----------
    file_features:
        the names of the bed file you want to load.
    chromosomes:
        chromosomes corresponding to the bed file. 
        If not specified, it's human by default
    path:
        if you want to specify the path where your bed file is. 
    sort: 
        if True, the bed file is sorted based on starting coordinates.
    """
    features_chrom = {}
    for c in chromosomes:
        features_chrom[c] = []
    with open(path + file_features) as features:
        for line in features:
            line = line.split("\t")
            if line[0][3:] in chromosomes:
                features_chrom[line[0][3:]].append([int(line[1]), int(line[2]), line[3]])
                
#    for c in chromosomes:
#        if features_chrom[c] ==  []:
#            del features_chrom[c]
#            
    if sort == True:
        for c in chromosomes:
            sorted(features_chrom[c], key=lambda x: x[0])
            
    return(features_chrom)

def make_windows(size, chromosomes=HUMAN, max_length=1000000000):
    """
    Generate windows/bins of the given size for the appropriate genome (default
    choice is human). 
    
    Parameters
    ----------
    size:
        size of the window you want
    chromosomes:
        Chromosomes of the species you are analysing
    max_length:
        maximum length given for every chromosome
    """
    features_chrom = {}
    start = range(0, max_length - size, size)
    end = range(size, max_length, size)
    for c in chromosomes:
        features_chrom[c] = [[start[i], end[i], ''.join(["chr_", c, "_", str(i)])] for i in range(len(end))]
    return(features_chrom)

###############################################################################
################################ PART 2 #######################################
###############################################################################

def extract_CG(input_sample, output_variable, head):
    """
    Extract relevant information for cytosines covered in a CG context
    It supposes that it is on the same format at the Methylpy/Ecker data
    Parameters
    ----------
    input_sample: file opened that contain the cytosines
    output_variable: the datatype that gather the methylation levels for a
        given type of features
    head: annotation of the file that we are reading (line to skip)
    """
    for cyt_called in input_sample:
        cyt_called = cyt_called.split("\t")
        if (cyt_called != head) and (cyt_called[3][0:2] == 'CG'):
            output_variable[cyt_called[0]].append([int(cyt_called[1]),
                       int(cyt_called[-3]), int(cyt_called[-2])])
    return(output_variable)
    
def extract_CH(sample, output_variable, head):
    """
    Extract relevant information for cytosines covered in a CH context
    It supposes that it is on the same format at the Methylpy/Ecker data
    Parameters
    ----------
    input_sample: file opened that contain the cytosines
    output_variable: the datatype that gather the methylation levels for a
        given type of features
    head: annotation of the file that we are reading (line to skip)
    """
    for cyt_called in sample:
        cyt_called = cyt_called.split("\t")
        if cyt_called != head and cyt_called[3][0:2] != 'CG':
            output_variable[cyt_called[0]].append([int(cyt_called[1]),
                       int(cyt_called[-3]), int(cyt_called[-2])])
    return(output_variable)


###############################################################################
    
#def extract_methylation (sample_name, feature, meth_type='CG', head=HEAD_Ecker, path='', chromosome=HUMAN):
#    
#    # It correspond to the annotation in the methylation call for methylpy.
#    # (at least for the Ecker dataset)
#    # I specify it here so I can automatically skip the line when I am reading the
#    # methylation files
#    
#    
#    # open the sample file
#    with open(path+sample_name) as sample:
#        
#        # initialize the variables
#        reduced_cyt = {} # to store cyt for each chrom (intermed output)
#        sample_bins = {} # the bins for each chrom (final output)
#        for c in chromosome:
#            reduced_cyt[c] = []
#            sample_bins[c] = [[0,0] for x in range(len(feature[c]))]
#
#        # specify the cytosine context you want to consider:
#        if meth_type == 'CG':
#            extract_CG(sample, reduced_cyt, head)
#        elif meth_type == 'CH':
#            extract_CH(sample, reduced_cyt, head)
#        elif meth_type == '':
#            for line in sample:
#                line = line.split("\t")
#                if line != head:
#                    reduced_cyt[line[0]].append([int(line[1]),int(line[-3]), int(line[-2])])
#        else:
#            print ("error in the choice of methylation context")
#            return()
#        return(reduced_cyt)
#        # so a this point, the cytosines that are relevant for me are stored 
#        # in the dict reduced_cyt
#        
#        # so one chromosome at a time, it would be nice to get the bins done. 
#        for c in chromosome:
#            cytosines = reduced_cyt[c] 
#            bins = feature[c]
#            curr_chrom = sample_bins[c]
#            i = 0
#            for j in range(len(bins)): # for every bins in a given chrom
#                # I am skipping cytosine that are before the beginning 
#                # of my current bin. 
#                while (i < len(cytosines)) and cytosines[i][0] < bins[j][0]:
#                    i += 1
#                # Once I got one cytosine that is after the beginning of my feature 
#                # I need to check if this feature is within the enhancer limits
#                # so if the position of the cytosine is not > to the end of the feature
#                if i<len(cytosines) and cytosines[i][0] <= bins[j][1]:
#                    curr_chrom[j][0] += cytosines[i][-2] # meth cyt
#                    curr_chrom[j][1] += cytosines[i][-1] # tot cyt
#
#                # check if the next cytosine fall into the current feature is important
#                # to this end, I have another pointer/iterator k. 
#                # at the next feature I will have to check from i but for the current 
#                # feature I need to check the next cytosine and I use the variable k for 
#                # this.
#                k = i+1
#                while k < len(cytosines) and cytosines[k][0] <= bins[j][1]:
#                    curr_chrom[j][0] += cytosines[k][-2] # meth cyt
#                    curr_chrom[j][1] += cytosines[k][-1] # tot cyt
#                    k += 1     
#       
#    return(sample_bins)
#    
###############################################################################  
###############################################################################  
    
def read_methylation_file(sample_name, meth_type, head, path, chromosome):
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

def methylation_level(reduced_cyt, feature, chromosome, threshold=1):
    """
    Measure the methylation for the feature you give as input using the reduce
    representation of the sample cytosine (output of read_methylation_file)
    Parameters
    ----------
    reduced_cyt: datatype that contained processed sample file. It only contains
        the cytosines that were in the genomic context you wanted to filter for.
        (output of read_methylation_file function).
    feature: the feature in the right datatype for which you want to determine
        the methylation level.
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes).
    """
    ## actually, to write sparse matrix I need a list, not a dictionary
    #meth_levels_bins = {key:[] for key in chromosome}
    meth_levels_bins = []

    for c in chromosome:
        meth_reads = np.zeros(len(feature[c]))
        tot_reads = np.zeros(len(feature[c]))
        nb_cyt = np.zeros(len(feature[c]))
        cytosines = reduced_cyt[c] 
        i = 0
        for j in range(len(feature[c])): # for every bins in a given chrom
            meth_reads = 0.0
            tot_reads = 1.0
            nb_cyt = 0
            # I am skipping cytosine that are before the beginning 
            # of my current bin. 
            while (i < len(cytosines)) and cytosines[i][0] < feature[c][j][0]:
                i += 1
            # Once I got one cytosine that is after the beginning of my feature 
            # I need to check if this feature is within the enhancer limits
            # so if the position of the cytosine is not > to the end of the feature
            if i<len(cytosines) and cytosines[i][0] <= feature[c][j][1]:
                meth_reads += cytosines[i][-2] # meth cyt read
                tot_reads += cytosines[i][-1] # tot cyt read
                nb_cyt += 1 # nb of cyt

            # check if the next cytosine fall into the current feature is important
            # to this end, I have another pointer/iterator k. 
            # at the next feature I will have to check from i but for the current 
            # feature I need to check the next cytosine and I use the variable k for 
            # this.
            k = i+1
            while k < len(cytosines) and cytosines[k][0] <= feature[c][j][1]:
                meth_reads += cytosines[i][-2] # meth cyt read
                tot_reads += cytosines[i][-1] # tot cyt read
                nb_cyt += 1  # nb of cyt
                k += 1
            ## actually, to write sparse matrix I need a list, not a dictionary
            if nb_cyt >= threshold:
                #meth_levels_bins[c].append(format(meth_reads/tot_reads, '.3f'))
                meth_levels_bins.append(format(meth_reads/tot_reads, '.3f'))
            else:
                #meth_levels_bins[c].append(np.nan)
                meth_levels_bins.append(np.nan)


    return(meth_levels_bins)

def write_not_sparse_meth(meth_to_write, name_file, writing_option='a', feature_names=None, cell_names=None):
    """
    takes the output of prep_methlevels and write everything in a given file.
    
    Parameters
    ----------
    """
    with open(name_file, writing_option) as file_to_write:
        if writing_option == 'w':
            file_to_write.write('%sparse matrix of methylation levels +1\n')
            file_to_write.write('%missing values are NaN not zeros.\n')
            if cell_names != None:
                file_to_write.write('%cell_names\t')
                file_to_write.write('\t'.join(cell_names))
                file_to_write.write('\n')
            if feature_names != None:
                file_to_write.write("feature_names\t")
                file_to_write.write('\t'.join(feature_names))
                file_to_write.write('\n')
    
        file_to_write.write('\t'.join([str(value) for value in meth_to_write]))
        file_to_write.write('\n')
        

### should I have an error if the size isn't right ?
    
    
def extract_methylation(sample_name, feature, meth_type=None, path='', #head=HEAD_Ecker,
	chromosome=HUMAN, threshold=1, col=None, write=False, writing_option='a', cell_names=None):
    """
    This function read the sample to produce a reduced representation of the
    sample (including only cytosines in a certain context) and then get
    corresponding methylation level in the feature given as input.
    Parameters
    ----------
    sample_name:
        name of the file to read to extract key information.
    feature:
    meth_type:
        CG, CH or not specified
    head: if there is header that you don't want to read. An annotation in the
        file you are reading. The default value is the Methylpy/Ecker header
    path: path of the to access the file of the sample you want to read. 
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes)
    cell_names is the list of cells you want to put as annotation at the beginning of your file
    """
    # It correspond to the annotation in the methylation call for methylpy.
    # (at least for the Ecker dataset)
    # I specify it here so I can automatically skip the line when I am reading the
    # methylation files
    ### doesn't work yet. You must have met reads and tot reads. 
    ## col argument contain the columns. We expect a list with the columns like such:
    # col = [chrom, pos_cyt, met_reads, unmet_reads, tot_reads]
    #If there is a missing column (either met_reads, unmer_reads, tot_reads) put 'NA'
    # also weird header ??
    if col != None:
    	chrom=col[0]
    	pos=col[1]
    	tot = col[3]
    	status=col[4]
    else:
    	chrom=0
    	pos=1
    	met=4
    	tot=5
    	status=3


    if meth_type == None:
    	reduced_cyt = read_meth_file(sample_name, path, chromosome, pos, met, tot, status)
    elif meth_type == 'CG':
    	reduced_cyt = read_meth_fileCG(sample_name, path, chromosome, pos, met, tot, status)
    elif meth_type == 'CH':
    	reduced_cyt = read_meth_fileCH(sample_name, path, chromosome, pos, met, tot, status)
    else:
    	# print a warning saying that the argument is not valid. We take all cytosines
    	reduced_cyt = read_meth_file(sample_name, head, path, chromosome)

    final_output = methylation_level(reduced_cyt, feature, chromosome, threshold)
    if write:
    	# I need to extract feature_names:
    	feature_names = extract_feature_names(feature)
    	write_not_sparse_meth(final_output, sample_name, writing_option, feature_names, cell_names)
    return(final_output)

## When I am making the feature names while generating features I should add the chromosome info.
def extract_feature_names(feature):
    bin_names = []
    for c in sorted(feature.keys()):
        bin_names += [x[-1] for x in feature[c]]

    return(bin_names)

def read_meth_fileCG(sample_name, path, chromosome, chrom=0, pos=1, met=4, tot=5, status=3):
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
    reduced_cyt = {key: [] for key in chromosome} # to store cyt for each chrom (intermed output)
    with open(path+sample_name) as sample:
        for line in sample:
            line = line.split('\t')
            if (line[status] in ['CGG', 'CGC', 'CGA', 'CGT']):
                reduced_cyt[line[chrom]].append((int(line[pos]), int(line[met]), int(line[tot])))
    return(reduced_cyt)

def read_meth_fileCH(sample_name, path, chromosome, chrom=0, pos=1, met=4, unmet=5, status=3):
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
    reduced_cyt = {key: [] for key in chromosome} # to store cyt for each chrom (intermed output)
    with open(path+sample_name) as sample:
        for line in sample:
            line = line.split('\t')
            if (line[status] not in ['CGG', 'CGC', 'CGA', 'CGT', 'mc_class']):
                reduced_cyt[line[chrom]].append((int(line[pos]), int(line[met]), int(line[unmet])))
    return(reduced_cyt)

def read_meth_file(sample_name, path, chromosome, chrom=0, pos=1, met=4, unmet=5, status=3):
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
    reduced_cyt = {key: [] for key in chromosome} # to store cyt for each chrom (intermed output)
    with open(path+sample_name) as sample:
        for line in sample:
            line = line.split('\t')
            if line[status] != 'mc_class':
                reduced_cyt[line[chrom]].append((int(line[pos]), int(line[met]), int(line[unmet])))
    return(reduced_cyt)
     
def filter_and_average_features(methylation_feature, feature, cell_name = '', split=False, threshold=1, average=''):
    """
    This function can do 3 different things :
        1- It is possible to filter features based on the number of cytosine covered per unit. 
        2- And/or to split the different feature contained. For example, a bed file that contain all
        TFBS can be consider as one type of feature or it can be split by the different TF.
        3- it can average the methylation level as a general average over all unit covered
        or it can average the average methylation level at each unit covered.
        
    Parameters
    ----------
    methylation_feature:    the data type that contain methylation level for the
        feature considered
    feature:    the feature data type for which you have the methylation level
    split:      do you want to split based on the feature names
    threshold:  the minimum number of cytosine per feature to keep it 
    average:    the type of average you want ('general', 'average', anything else mean
    no average done)
    """
    final_output = [cell_name]
    list_methylation = make_list(methylation_feature)
    if split==True:
        list_feature = [x[-1] for x in make_list(feature)] 
    else:
        list_feature = ['standard_name' for x in range(len(list_methylation))]
    name_feature = list(set(list_feature))
    output_average = {}
    
    if average == 'general':
        for keys in name_feature:
            output_average[keys] = [0,0]
        for index in range(len(list_methylation)):
            if list_methylation[index][-1] >= threshold:
                output_average[list_feature[index]][0] += list_methylation[index][0]
                output_average[list_feature[index]][1] += list_methylation[index][1]

    elif average == 'average':
        for keys in name_feature:
            output_average[keys] = []
        for index in range(len(list_methylation)):
            if list_methylation[index][-1] >= threshold and list_methylation[index][1]>0:
                output_average[list_feature[index]].append(list_methylation[index][0]/list_methylation[index][1])
        #return(output_average)
    # can be optimized by doing 2 arrays and then divide every elements
                
    else:
        for keys in name_feature:
            output_average[keys] = []
        for index in range(len(list_methylation)):
            if list_methylation[index][-1] >= threshold:
                output_average[list_feature[index]].append(list_methylation[index][0:2])
        return(output_average)
        
    for key in sorted(output_average.keys()):
        if average == 'average':
            final_output.append(np.average(output_average[key]))
        elif average == 'general' and output_average[key][1]>0:
            final_output.append(float(output_average[key][0])/output_average[key][1])
        else:
            final_output.append("nan")
        
    return(final_output)
    
def filter_and_average_features_chrm(methylation_feature, feature, cell_name = '', split=False, threshold=1, average=''):
    """
    This function can do 3 different things :
        1- It is possible to filter features based on the number of cytosine covered per unit. 
        2- And/or to split the different feature contained. For example, a bed file that contain all
        TFBS can be consider as one type of feature or it can be split by the different TF.
        3- it can average the methylation level as a general average over all unit covered
        or it can average the average methylation level at each unit covered.   
    Parameters
    ----------
    methylation_feature:    the data type that contain methylation level for the
        feature considered
    feature:    the feature data type for which you have the methylation level
    split:      do you want to split based on the feature names
    threshold:  the minimum number of cytosine per feature to keep it 
    average:    the type of average you want ('general', 'average', anything else mean
    no average done)
    """
    final_output = [cell_name]
    for chrm in methylation_feature.keys():
        list_methylation = [x for x in methylation_feature[chrm]]
        if split==True:
        	list_feature = [x[-1] for x in make_list(feature)] 
        else:
        	list_feature = ['standard_name' for x in range(len(list_methylation))]
        name_feature = list(set(list_feature))
        output_average = {}
        if average == 'general':
            for keys in name_feature:
                output_average[keys] = [0,0]
            for index in range(len(list_methylation)):
                if list_methylation[index][-1] >= threshold:
                    output_average[list_feature[index]][0] += list_methylation[index][0]
                    output_average[list_feature[index]][1] += list_methylation[index][1]

        elif average == 'average':
            for keys in name_feature:
                output_average[keys] = []
            for index in range(len(list_methylation)):
                if list_methylation[index][-1] >= threshold  and list_methylation[index][1]>0:
                    output_average[list_feature[index]].append(list_methylation[index][0]/list_methylation[index][1])
        #return(output_average)
    # can be optimized by doing 2 arrays and then divide every elements
        else:
            for keys in name_feature:
                output_average[keys] = []
            for index in range(len(list_methylation)):
                if list_methylation[index][-1] >= threshold:
                    output_average[list_feature[index]].append(list_methylation[index][0:2])
            return(output_average)

        for key in sorted(output_average.keys()):
            if average == 'average':
                final_output.append(np.average(output_average[key]))
            elif average == 'general' and output_average[key][1]>0:
                final_output.append(float(output_average[key][0])/output_average[key][1])
            else:
                final_output.append("nan")

    return(final_output)

def prep_methlevels(methylation_feature, name_cell, threshold=1):
    """
    Takes methylation levels, replace bins that are not covered enough by NaN,
    otherwise get the average methylation level for the bin and return the final
    values in alist. Ready to write. 
    
    Parameters
    ----------
    """
    output_list = [name_cell]
    for element in make_list(methylation_feature):
        if element[2] >= threshold:
            output_list.append(float(element[0])/element[1])
        else:
            output_list.append(NaN)
    return(output_list)
  
def write_list(the_list, file_to_write):
    for value in the_list:
        file_to_write.write(str(value))
        file_to_write.write("\t")
    file_to_write.write("\n")
    return()
    
def write_methlevel(meth_to_write, name_file, cell, writing_option='a', feature_names=None):
    """
    takes the output of prep_methlevels and write everything in a given file.
    
    Parameters
    ----------
    """
    file_to_write = open(name_file, writing_option)

    if feature_names != None:
        feature_to_write = ["sample_name"]
        feature_to_write += sorted(feature_names)
        write_list(feature_to_write, file_to_write) # should add an condition
                                            # to be sure features is a dictionary
    meth_to_write = [cell] + meth_to_write
    write_list(meth_to_write, file_to_write)
    file_to_write.close()
    return()

###############################################################################
################################### PART 3 ####################################
###############################################################################

def summary_stat1(all_data_features, sample_name, nb ,chromosomes = HUMAN, write=False, f=''):
    """
    old function (not necessarily relevant anymore)
    get the coverage of a sample for the feature covered.
    """
    covered = 0
    tot = 0
    for c in HUMAN:
        for element in all_data_features[sample_name][c]:
            tot += 1
            if element[1] != 0:
                covered +=1

    print(nb, '  -  ', sample_name, '  -  ', covered*100/float(tot))
    if write == True:
        t = ''.join([nb, '  -  ', sample_name, '  -  ', covered*100/float(tot)])
        #t = str(t)
        f.write(t)
    #return()    
###############################################################################
    
# some statistics about the TFBS coverage and methylation level
def summary_stat2(all_data_features, chromosomes = HUMAN, write=False, f=''):
    """
    old function (not necessarily relevant anymore)
    for TFBS summary statistics. Output meth level per TFBS (if not already an average) 
    nb of cytosines per TFBS ??
    """
    tot_cyt_TFBS = {}
    meth_level = {}
    meth_cyt_count = {}
    for sample, chroms in all_data_features.items():
        print(sample[0:-1])
        meth_level[sample] = []
        meth_cyt_count[sample] = []
        tot_cyt_TFBS[sample] = []
        for c in HUMAN:
            for TFBS_bin in chroms[c]:
                if TFBS_bin[1] != 0:
                    meth_cyt_count[sample].append(TFBS_bin[0])
                    tot_cyt_TFBS[sample].append(TFBS_bin[1])
                    meth_level[sample].append(float(TFBS_bin[0])/TFBS_bin[1])
    
        print('meth level per TFBS ', np.mean(meth_level[sample]))
        print('nb of cyt per TFBS ', np.mean(tot_cyt_TFBS[sample]))
        print('\n')
        if write == True:
            f.write(np.mean(meth_level[sample]))
            f.write(np.mean(tot_cyt_TFBS[sample]))
            f.write('\n')
    
    return(meth_level)
    
###############################################################################    
    
def plot_meth_level(meth_level, title='', label_x='', label_y='', pdf_file = '', shape=10, k=1):
    """
    old function (not necessarily relevant anymore)
    plot methylation level of samples for different features. 
    """
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    fig.set_size_inches(18.5, 10.5)
  
    j = 1
    for sample, values in meth_level.items():
        title2 = sample+ title
        plt.subplot(k, shape, j)
        sns.distplot(values)
        plt.title(title)
        plt.xlabel(label_x)
        plt.ylabel(label_y) 
        fig.suptitle(title2)
        j += 1
    
    plt.show()
    if pdf_file != '':
        #plt.savefig(pdf_file, format='pdf')
        pdf_file.savefig(fig)

#    return(fig)
###############################################################################
################################# PART 4 ######################################
###############################################################################
def read_matrix(matrix_file):
    """ 
    just a function to read a file containing a matrix that we want to use
    to cluster and plot etc.
    BONUS : I will return an AnnData like matrix. So it is directly usable for
    scanpy after that.
    """
    raw_matrix = []
    matrix = [] # np.ndarray
    variables = [] # pd.DataFrame, dict --> one dim. of length #variables
    obs = [] # pd.DataFrame, dict --> one dim. of length #observations
    with open(matrix_file) as m:
        raw_matrix = m.readlines()
    variables = raw_matrix[0].split('\t')[1:]
    for sample in raw_matrix[1:]:
        line = sample.split('\t')
        obs.append(line[0])
        matrix.append([float(x) for x in line[1:-1]])
    return matrix, variables, obs
   
###############################################################################
    
def make_list(dico):
    """
    Make information stored in  dictionary with chromosomes as key
    into one list. (lost of chromosomal information)
    """
    final_list = []
    for c in sorted(dico.keys()):
        final_list += dico[c]
    return final_list

def mean_meth(list_bins):
    """
    take the list containing the feature units and calculate the methylation
    level for the bins that are covered and return the average methylation
    level over the different covered units.
    """
    meth_level = []
    for unit in list_bins:
        if unit != [0, 0, 0]:
            meth_level.append(float(unit[0])/unit[1])
    return np.mean(meth_level)

###############################################################################

def readandimputematrix(file_name, min_coverage=1):
    with open(file_name) as f:
        file = f.readlines()

    # separate annotation from data    
    head_var = file[0]
    head_var = head_var.split('\t')
    # Then, extract the sample names
    sample_names = []
    data_raw = []
    for l in file[1:]:
        l = l.split('\t')
        sample_names.append(l[0])
        data_raw.append(l[1:])

    # clear memory of useless variables 
    del file
    
    ##########################################
    # now, removing empty columns
    empties = []
    partial = []
    full =  []
    for index in range(1, len(data_raw[0])):
        column = [element[index] for element in data_raw]
        if len(list(set(column))) == 1:
            empties.append(index)
        elif len(list(set(column))) <= min_coverage:
            partial.append(index)
        else:
            full.append(index)
         
    ##########################################
    intermed_matrix = []
    name_windows_covered = []
    # let's remove the compltetly uninformative columns
    for index in range(1, len(head_var[1:])):
        if index in full:
            intermed_matrix.append([element[index] for element in data_raw])
            name_windows_covered.append(head_var[index])

    ########################################
    # imputing values.
    imputed_matrix = []
    for row in intermed_matrix:
        imputed_row = []
        if "nan" in row:
            mean = np.mean([float(e)  for e in row if e != "nan"])
            for element in row:
                if element == "nan":
                    imputed_row.append(str(mean))
                else: 
                    imputed_row.append(element)
            imputed_matrix.append(imputed_row)
        else:
            imputed_matrix.append(row)
        
    imputed_matrix = np.matrix(imputed_matrix).transpose()
    return(imputed_matrix, sample_names, name_windows_covered)



def alternative_windows(annotation_file, path_file="",  genome=MOUSE, minimum_annot=2):
    """
    It creates windows delimited by the annotation value (single base pair considered)
  
    Parameters
    ----------
    annotation_file:    files containing annotation coordinates to delimit the
                        windows. bed files (or similar).
    path_file:    path to the annotation file
    genome:      genome to consider
    minimum_annot:    if you want to filter out chromosomes that don't have enough
                      annotations to generate windows. 
    """
    with open(path_file+annotation_file) as f:
        raw_annot = f.readlines()
    raw_annot = [x.split('\t') for x in raw_annot]
    
    annot_windows = {key: [] for key in MOUSE}
    for line in raw_annot:
        annot_windows[line[0][3:]].append(int(line[1]))
     
    key_to_del = []
    for key, value in annot_windows.items():
        if len(value) < minimum_annot:
            key_to_del.append(key)
    if key_to_del != []:
        print("deleted chromosomes:")
        for key in key_to_del:
            print("chr", key, len(annot_windows[key]))
            del annot_windows[key]

    new_annot_windows = {key : [] for key in annot_windows.keys()}
    for key in sorted(annot_windows.keys()):
        tmp_window = [0]
        ## along the chromosome
        for v in annot_windows[key]:
            tmp_window.append(v)
            tmp_window.append("_".join(["chr", key, str(tmp_window[0]), str(v)]))
            new_annot_windows[key].append(tmp_window)
            tmp_window = [v]

        # to end a chromosome
        tmp_window.append(float("inf"))
        tmp_window.append("_".join(["chr", key, str(tmp_window[0]), "inf"]))
        new_annot_windows[key].append(tmp_window)
        print(tmp_window)
        tmp_window = [0]
    return(new_annot_windows)

def names_annot(loaded_annot):
    names = []
    for key in sorted(loaded_annot.keys()):
        for v in loaded_annot[key]:
            names.append(v[-1])
    return(names)
