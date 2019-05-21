import anndata as ad

def load_metadata(adata, metadata_file, path=''):
    """
    Load observational metadata in adata.obs.
    Input metadata file as csv and adata object.

    first raw of the metadata file is considered as a header

    Paramters
    ---------

    adata: initial AnnData object

    metadata_file: csv file containing as a first column the cell names and in the
        rest of the columns any king of metadata to load
    
    Return
    ------
    Annotated AnnData
    """
    dict_annot = {}
    with open(path+metadata_file) as f:
        head = f.readline().split(';')
        file = f.readlines()
    for key in head:
        dict_annot[key] = []
    data = [line.split(';') for line in file]

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
