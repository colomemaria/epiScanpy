# chromosomes for 2 principal species. If you work with another genome
# the chromosomes will have to be specified
HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2',
         '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2',
         '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']
# it assume that you take a bed file as impute
# it save the features position and name like [start, end, name]
# in a dictionary that has chromosomes as keys. 
def load_features(file_features, chromosomes=HUMAN, path="", sort=False):
    """
    Description function. 
    
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
    Description function. 
    
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
