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
