import anndata as ad

def load_metadata(adata, metadata_file, path='', separator=';'):
    """
    Load observational metadata in adata.obs.
    Input metadata file as csv/txt and the adata object to annotate.

    first raw of the metadata file is considered as a header
    first column contain the cell name

    Paramters
    ---------

    adata: initial AnnData object

    metadata_file: csv file containing as a first column the cell names and in the
        rest of the columns any king of metadata to load

    path: pathe to the metadata file

    separator: ';' or "\t", character used to split the columns
    
    Return
    ------
    Annotated AnnData 
    """
    dict_annot = {}
    with open(path+metadata_file) as f:
        head = f.readline().split(separator)
        file = f.readlines()
    for key in head:
        dict_annot[key] = []
    data = [line.split(separator) for line in file]

    for name in adata.obs_names:
        name = name.split('.')[0]
        found = False
        for line in data:
            if name == line[0]:
                i = 0
                for key in head:
                    dict_annot[key].append(line[i])
                    i += 1 
                found = True
                continue
        # if we could not find annotations
        if found == False:
            for key in head:
                dict_annot[key].append('NA')

    for key in head:
        adata.obs[key] = dict_annot[key]

    return()