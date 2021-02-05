import anndata as ad
import pandas as pd

def load_metadata(adata, metadata_file, path='', separator=';', remove_index_str = None):
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
    
    remove_index_str: a list of string to be removed in the index of AnnData object.
    Default value is None. For example, if the index is ['/path/to/file1.txt','/path/to/file2.txt']
    and remove_index_str = ['/path/to/','.txt'], then the index of AnnData object 
    will be changed to ['file1','file2']
    
    Return
    ------
    Annotated AnnData 
    """
    # dict_annot = {}
    # with open(path+metadata_file) as f:
    #     head = f.readline().strip().split(separator)
    #     file = f.readlines()
    # for key in head:
    #     dict_annot[key] = []
    # data = [line.strip().split(separator) for line in file]
    # for name in adata.obs_names:
    #     # this line is not always true. It depends on how the format of the data are
    #     name = name.split('.')[0]
    #     found = False
    #     for line in data:
    #         if name == line[0]:
    #             i = 0
    #             for key in head:
    #                 dict_annot[key].append(line[i])
    #                 i += 1 
    #             found = True
    #             continue
    #     # if we could not find annotations
    #     if found == False:
    #         for key in head:
    #             dict_annot[key].append('NA')
    # for key in head:
    #     adata.obs[key] = dict_annot[key]
    metadata = pd.read_csv(path+metadata_file, sep = "\t", header = 0)
    str_index = adata.obs.index
    if remove_index_str:
        for value in remove_index_str:
            str_index = str_index.str.replace(value,'',regex=False)
    df = pd.DataFrame('NA', index=str_index, columns=metadata.columns)
    for key,value in metadata.iterrows():
        try:
            df[value.index[0]][value[0]] =  value[0]
            df[value.index[1]][value[0]] =  value[1]
        except:
            continue
    adata.obs.index = str_index
    for key in df.columns:
        adata.obs[key] = df[key]
