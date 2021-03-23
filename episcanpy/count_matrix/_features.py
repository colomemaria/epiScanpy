import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ..preprocessing._episcanpy_mo_fcts import load_gtf_file

# chromosomes for 2 principal species. If you work with another genome
# the chromosomes will have to be specified
# mitochondrial genome not included
MOUSE = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19','X', 'Y']
MOUSE_SIZE = [195200000, 181800000, 159800000, 156900000, 151800000, 149650000, 145050000, 130200000, 124400000, 130600000,  
        122000000, 120150000, 120950000, 125250000, 104150000, 98050000, 95350000, 90800000, 61500000, 169550000, 91500000]

HUMAN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
        '21', '22','X', 'Y']        
HUMAN_SIZE = [249000000, 242250000, 198350999, 190250000, 181600000, 170850000, 159400000, 145200000, 138450000, 133850000, 
        135150000, 133350000, 114400000, 107100000, 102050000, 90400000, 83300000, 80400000, 58650000, 64500000,
        46750000, 50850000, 156100000, 57300000]


def filter_df(df, column_name, filter_criteria_list):
    filter_annot = []
    for value in df[column_name]:
        if value in filter_criteria_list:
            filter_annot.append('keep')
        else:
            filter_annot.append('discard')
    df['filter'] = filter_annot  
    df = df[df['filter']=='keep']
    del df['filter']
    return(df)


def load_features_bed(file_features, chromosomes=HUMAN, path=""):
    """
    load features when input file is bed
    """
    features_chrom = {}
    for c in chromosomes:
        features_chrom[c] = []
    with open(path + file_features) as features:
        for line in features:
            line = line.strip().split("\t")
            if line[0][3:] in chromosomes:
                try:
                    features_chrom[line[0][3:]].append([int(line[1]), int(line[2]), line[3]])
                except:
                    features_chrom[line[0][3:]].append([int(line[1]), int(line[2]) ])
                         
    return(features_chrom)

def load_features_gff(file_features,
                      chromosomes=None, # can be str or a list of str
                      filter_per_source=None, # can be str or a list of str
                      filter_per_feature_type=None, # can be str or a list of str
                      path="",
                     sort=False):
    """
    load features when input file is gff. 
    
    Using chromosomes, filter_per_source and/or filter_per_feature_type it is
    possible to load only a subset of the features. 
    
    Parameters
    ----------
    file_features:
        the names of the bed file you want to load.
    chromosomes:
        chromosomes corresponding to the bed file. 
        If not specified, it's human by default
    path:
        if you want to specify the path where your bed file is. 
        
    input_file_format:
        if None, the input format is deduced from the extension name (admitted bed, gtf, gff)
        if str specified, it overlook the rxtension name and load the feature file as the 
        specified input file format.
    sort: 
        if True, the bed file is sorted based on starting coordinates.
    
    """
    # load the gff file
    df = []
    with open(path+file_features) as f:
        for line in f: 
            if line[0] != '#':
                df.append(line.rstrip('\n').split('\t'))
                

    columns = ['sequence', 'source', 'feature', 'start', 'end',
               'score', 'strand', 'phase', 'attribute']
    df = pd.DataFrame(np.array(df), columns=columns)
    
    # filter the df
    if chromosomes != None:
        if type(chromosomes) == str:
            chromosomes = [chromosomes]
        df = filter_df(df,
                       column_name='sequence',
                       filter_criteria_list=chromosomes)
        
    if filter_per_source != None:
        if type(filter_per_source) == str:
            filter_per_source = [filter_per_source]
        df = filter_df(df,
                       column_name='source',
                       filter_criteria_list=filter_per_source)
        
    if filter_per_feature_type != None:
        if type(filter_per_feature_type) == str:
            filter_per_feature_type = [filter_per_feature_type]
        df = filter_df(df,
                       column_name='feature',
                       filter_criteria_list=filter_per_feature_type)
        
    # convert the df into the right format to build a count matrix from this feature set
    feature_output = {key: [] for key in set(df['sequence'])}
    index = 0
    for line in df.values:
        feature_output[line[0]].append([int(line[3]),
                                        int(line[4]),
                                        "_".join([line[2], line[1], str(index)])
                                       ])
        index += 1
    del df
        
    if sort == True:
        for c in chromosomes:
            sorted(features_chrom[c], key=lambda x: x[0])
            
    return(feature_output)

def load_features_gtf(file_features,
                      chromosomes=None, # can be str or a list of str
                      filter_per_source=None, # can be str or a list of str
                      filter_per_feature_type=None, # can be str or a list of str
                      path="",
                      sort=False):
    """
    oad features when input file is gtf. 
    
    Using chromosomes, filter_per_source and/or filter_per_feature_type it is
    possible to load only a subset of the features. 
    
    """
    df = load_gtf_file(path+file_features)
    
    # filter the df
    if chromosomes != None:
        if type(chromosomes) == str:
            chromosomes = [chromosomes]
        df = filter_df(df,
                       column_name='seqname',
                       filter_criteria_list=chromosomes)
        
    if filter_per_source != None:
        if type(filter_per_source) == str:
            filter_per_source = [filter_per_source]
        df = filter_df(df,
                       column_name='source',
                       filter_criteria_list=filter_per_source)
        
    if filter_per_feature_type != None:
        if type(filter_per_feature_type) == str:
            filter_per_feature_type = [filter_per_feature_type]
        df = filter_df(df,
                       column_name='feature',
                       filter_criteria_list=filter_per_feature_type)
        
    # convert the df into the right format to build a count matrix from this feature set
    feature_output = {key: [] for key in set(df['seqname'])}
    index = 0
    for line in df.values:
        feature_output[line[0]].append([int(line[3]),
                                        int(line[4]),
                                        "_".join([line[2], line[1], str(index)])
                                       ])
        index += 1
        
    del df
    
    if sort == True:
        for c in chromosomes:
            sorted(features_chrom[c], key=lambda x: x[0])
            
    return(feature_output)
    

def load_features(file_features, chromosomes=HUMAN, path="", input_file_format=None, sort=False):
    """
    The function load features is here to transform a bed file into a usable 
    set of units to measure methylation levels. 
    It has to be a bed-like file. You also need to specify the chromosomes you 
    use as a list of characteres like ['1', '7', '8', '9', 'X', 'Y', 'M'].
    The chromosome list you give as input can be not ordered.
    If you don't specify the chromosomes, the default is the human genome
    (including X, Y and mitochondrial DNA).
    THE BED (and gtf) FILE need to be sorted. (upcoming an option to sort the
    file).
    
    The output is a dictionary where the keys are chromosomes and the value is
    a list containing [start, end, name] for every feature extracted.
    
    This function will load the entire annoation file. If you want to use only parts of gtf/gff files
    please look the functions load_features_gff and load_features_gtf.
    
    Parameters
    ----------
    file_features:
        the names of the bed file you want to load.
    chromosomes:
        chromosomes corresponding to the bed file. 
        If not specified, it's human by default
    path:
        if you want to specify the path where your bed file is. 
        
    input_file_format:
        if None, the input format is deduced from the extension name (admitted bed, gtf, gff)
        if str specified, it overlook the rxtension name and load the feature file as the 
        specified input file format.
    sort: 
        if True, the bed file is sorted based on starting coordinates.
    """
    
    if input_file_format==None:
        input_file_format = file_features[-3:]
             
    if input_file_format not in ['bed', 'gtf', 'gff']:
        # return warning that the input format isn't correct and that you need to either 
        # provide the correct format or specify the format with the input_file_format parameter
        print('warning')
    elif input_file_format=='bed':
        features_chrom = load_features_bed(file_features, chromosomes, path)
    elif input_file_format=='gtf':
        features_chrom = load_features_gtf(file_features, chromosomes, path)
    elif input_file_format=='gff':
        features_chrom = load_features_gff(file_features, chromosomes, path)
    else:
        print(input_file_format)
        
    if sort == True:
        for c in chromosomes:
            sorted(features_chrom[c], key=lambda x: x[0])
            
    return(features_chrom)

def make_windows(size, chromosomes = 'human', chromosome_sizes = None, max_length=1000000000):
    """
    Generate windows/bins of the given size for the appropriate genome (default
    choice is human). 
    
    Parameters
    ----------
    size:
        size of the window you want
    chromosomes:
        Chromosomes of the species you are analysing. Pre-defined chromosomes are 'human' and 'mouse'. Default value is 'human'.
        For other oranisms, chromosomes should be defined as, for example, chromosomes = ['1', '2', ... , 'X', 'Y'].
    chromosome_sizes: 
        an array for chromosome sizes. If chromosomes is set to 'human' or 'mouse', chromosome_sizes will be ignored. 
        For other oranisms, chromosome_sizes must have the same length as chromosomes and should be defined as, for example, 
        chromosomes = ['1', '2', ... , 'X', 'Y'].
    max_length:
        maximum length given for every chromosome. For other oranisms, when chromosomes is assigned and chromosome_sizes = None,
        max_length will be used instead.
    """
    features_chrom = {}
    
    if chromosomes == 'human':
        chromosomes = HUMAN
        chromosome_sizes = HUMAN_SIZE
    elif chromosomes == 'mouse':
        chromosomes = MOUSE
        chromosome_sizes = MOUSE_SIZE
    elif isinstance(chromosomes, list):
        if isinstance(chromosome_sizes, list):
            if len(chromosomes) != len(chromosome_sizes):
                print("Error: the lengths of chromosomes and chromosome_sizes are not matched.")
        else:
            chromosome_sizes = [max_length for i in range(len(chromosomes))]
    else:
        print("Error: chromosomes must be assigned as 'human', 'mouse', or a list")
        
    for c in range(len(chromosomes)):
        start = range(1, chromosome_sizes[c] - size, size)
        end = range(size, chromosome_sizes[c], size)
        features_chrom[chromosomes[c]] = [[start[i], end[i], ''.join(["chr", chromosomes[c], "_", str(start[i]), "_", str(end[i])])] for i in range(len(end))]
        
    return(features_chrom)


def size_feature_norm(loaded_feature, size):
    """If the features loaded are too smalls or of different sizes, 
    it is possible to normalise them to a unique given size by extending the 
    feature coordinate in both directions.

    Parameters
    ----------
    loaded_feature: loaded feature to normalise

    size: desired size of the feature

    Return
    ------
    Update the input features

    """
    for key in loaded_feature.keys():
        for i in range(len(loaded_feature[key])):
            length = loaded_feature[key][i][1] - loaded_feature[key][i][0]
            add_length = size - length
            if add_length%2 == 1:
                loaded_feature[key][i][0] = loaded_feature[key][i][0] - add_length//2 -1
                loaded_feature[key][i][1] = loaded_feature[key][i][1] + add_length//2
            else:
                loaded_feature[key][i][0] = loaded_feature[key][i][0] - add_length//2
                loaded_feature[key][i][1] = loaded_feature[key][i][1] + add_length//2
    
def plot_size_features(loaded_feature, bins=50, return_length=False):
    """
    Plot the different feature sizes in an histogram. 
    """
    length = []
    for key in loaded_feature.keys():
        for line in loaded_feature[key]:
            length.append(line[1] - line[0])
        
    np.histogram(length)
    plt.hist(length, bins)
    plt.show
    if return_length:
        return(length)
    
    
def name_features(loaded_features):
    """
    Extract the names of the loaded features, specifying the chromosome they originated from.
    It also contain the feature coordinates and an unique identifier.
    """
    feat_names = []
    i = 0
    for c in loaded_features.keys():
        for name in loaded_features[c]:
            #add_name = '_'.join(['chr', c, name[-1].rstrip('\n'), str(i)])
            add_name = '_'.join([c, str(name[0]), str(name[1])])
            add_name ='chr' + add_name
            # the feature names will be like chr1_1234_1245
            if add_name[-1] =='\n':
                add_name = add_name[:-1]
            feat_names.append(add_name)
            i += 1
    return(feat_names)
