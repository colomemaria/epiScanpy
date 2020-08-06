# enables the %%R magic
#%load_ext rpy2.ipython

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
readRDS = robjects.r['readRDS']

import pandas as pd
from warnings import warn


import episcanpy.api as epi
import scanpy as sc
import anndata as ad

import anndata2ri
anndata2ri.activate()

# import R's "SnapATAC" package
snapatac = importr('SnapATAC')

def digRS4(RS4, verbose=True):
    """
    Extract elements of an RS4 object.
    
    """
    output_list = [type(RS4), RS4.rclass, RS4.slots.keys()]
    
    if verbose:
        print('running function digRS4')
        print(type(RS4), RS4.rclass)
        print('elements of the RS4', RS4.slots.keys())
        
    # go through the different elements, extract the data and store it in a list
    for i in tuple(RS4.slots.keys()):
        list_type = []
        type_i = type(RS4.slots[i]) # extract the type of the element extracted. 
        if verbose:
            print(tuple(RS4.rclass), '\t', i, type(RS4.slots[i]))
        # if the element is not an RS4 you simply store it in the list
        if type_i != rpy2.robjects.methods.RS4:
            list_type.append(RS4.slots[i])
        # if the element is an RS4 we call this function to go one step deeper 
        ### WARNING: need to check it is not and endless loop
        else:
            if verbose:
                print('ERROR: RS4 object containing RS4')
            list_type.append(digRS4(RS4.slots[i], verbose=verbose))
        output_list.append(list_type)
    return(output_list)


def extract_snap(rds_file_name, verbose=True):
    """
    reads an rds file containing a snap object. 
    Outputs an anndata matrix
    
    Warning:
    the current version do not transfer Granges information
    
    Input:
    ------
    
    rds_file_name --> file name, it had to be an rds file
    
    verbosity --> If True, prints information regarding the different steps to read
    and convert thesnap object into an anndata.
    
    save --> if a file name specified, it will save the file as h5ad
    
    Output:
    -------
    Anndata object
    
    """
    #### Load the snao object and check what's in it
    if verbose:
        print('READ THE RDS FILE:\n')
    snap = readRDS(rds_file_name)
    
    
    if verbose:
        # print the type and the class of snap
        print(type(snap), tuple(snap.rclass))
        # print keys
        print('keys: \t', tuple(snap.slots.keys()))
    if 'snap' not in tuple(snap.rclass):
        warn('ERROR, the input file needs to be an rds file containing a snap object\n')
        return()
    elif verbose:
        print('continue with the snap object.\n\n\n')
    
    ############
    ###  extract the different elements of the snap object:
    if verbose:
        print('EXTRACT ELEMENTS COMPOSING THE SNAP OBJECT:\n')
        
    #extract all elements to export into the Anndata that are not RS4 themself
    elements = {}
    # go through the individual elements of the snap object
    for n in tuple(snap.slots.keys()):
        if verbose:
            print(n, type(snap.slots[n]))
            
        #######
        # the case of RS4 and GRanges
        if type(snap.slots[n]) == rpy2.robjects.methods.RS4 and ('GRanges' in snap.slots[n].rclass):
            
            tmp_list = [type(snap.slots[n]), tuple(snap.slots[n].rclass)]
            warn(''.join(['\n', n, 'partially excluded from the conversion.\n GRanges are not currently fully transfered into the Anndata.\n']))
            granges = snap.slots[n]
            tmp_list.append(tuple(granges.slotnames()))
            
            tmp_dict={}
            for m in tuple(granges.slotnames()):
                if verbose:
                    print(m, type(granges.slots[m]))
                # only salvage the elements that are not RS4
                if type(granges.slots[m]) != rpy2.robjects.methods.RS4:
                    tmp_dict[m] = granges.slots[m]
                else:
                    tmp_list2 = []
                    tmp_feat = granges.slots[m]
                    for p in tuple(tmp_feat.slotnames()):
                        if verbose:
                            print('\t', p, type(tmp_feat.slots[p]))
                        tmp_list2.append([p, type(tmp_feat.slots[p]), tmp_feat.slots[p]])
                    
                    tmp_dict[m] = [tuple(granges.slots[m].slotnames()), tmp_list2]
                tmp_list.append(tmp_dict)

            elements[n] = tmp_list
            del tmp_list, tmp_list2, tmp_dict
            if verbose:
                print('\n')
        
        #######
        # the case of RS4 (not GRanges)
        # if the object is an RS4 we will call a function toextract the subparts:
        elif type(snap.slots[n]) == rpy2.robjects.methods.RS4:
            if verbose:
                print('\nRS4 object: ', n, tuple(snap.slots[n].rclass), tuple(snap.slots[n].slotnames()), '\n')
                
            ### call function
            elements[n] = digRS4(RS4=snap.slots[n], verbose=verbose)
            
            
        #######
        # the rest of the cases
        # if it is not an RS4 we simpoly convert it to a python data type and store it in the dictionary
        else:
            elements[n] = [type(snap.slots[n]), snap.slots[n]]
            
            
    return(elements)


def make_Anndata(input_data,
                 mtx_name='bmat',
                 feature_names=None,
                 metadata_name='metaData', 
                 jaccard=True,
                 jaccard_key='jmat',
                 save_all_rds=False,
                 save=None,
                 copy=True):
    """
    Currently, having multiple layers isn't possible.
    It is possible to add the jaccard matrix.
    
    """
    
    adata = ad.AnnData(X=input_data[mtx_name][1],
            obs=input_data['metaData'][1])

    if feature_names != None:
        if 'GRanges' in input_data[feature_names][1]:
            adata.var = input_data[feature_names][3]['elementMetadata']
            
    ## extra annotations ?
    extra_annot = ['barcode', 'file', 'cluster', 'sample']
    for extra in extra_annot:

        if (extra in input_data.keys()) and (extra not in adata.obs.columns) and (len(input_data[extra][1].tolist()) == len(adata.obs_names.tolist())):
            adata.obs[extra] = input_data[extra][1].tolist()
            
        
    
    # save the class snap:
    if 'class' in input_data.keys():
        adata.uns['class'] = input_data['class'][1]
    
    ## Add tsne, umap
    for n in ['tsne', 'umap']:
        if n in input_data.keys():
            name = ''.join(['X_',n,'_snap'])
            if input_data[n][1].shape == (len(adata.obs_names.tolist()), 2):
                adata.obsm[name] = input_data[n][1]
            else:
                warn(''.join(['WARNING: the dimension of ',
                              n,
                              ' does not match the expected shape. Therefore it is stored in adata.uns\n']))
                adata.uns[name] = input_data[n]
                
            
    if jaccard and (jaccard_key in input_data.keys()):
        
        ## check that the jaccard distance was computed properly
        if mtx_name not in input_data[jaccard_key][8][0]:
            warn(''.join(['WARNING: the jaccard distance is not computed using ', mtx_name, '\n',
                 'The jaccard distance was computed for: ', str(input_data[jaccard_key][8][0]), '\n']))
            
        adata.uns['jaccard'] = {'index':input_data[jaccard_key][3][0],
                                'normalised':input_data[jaccard_key][4][0],
                                'params':{'mtx': input_data[jaccard_key][8][0],
                                          'vector1': input_data[jaccard_key][5][0],
                                          'vector2': input_data[jaccard_key][6][0]
                                         }
                               }
    # save the rest in uns
    if save_all_rds:
        for n in input_data.keys():
            adata.uns[n] = input_data[n]
         
    if save!=None:
        save = save.rstrip('.h5ad')
        adata.write('.'.join([save, 'h5ad']))
        
    if copy:
        return(adata)
    
    
    
def snap2anndata(rds_file_name,
                 path='',
                 mtx_name='bmat',
                 feature_names='feature',
                 metadata_name='metaData', 
                 jaccard=True,
                 jaccard_key='jmat',
                 save_all_rds=False,
                 verbose=False,
                 copy=True,
                 save=None):
    """
    Convert a snap object saved as .rds into an Anndata object.
    
    In the current version of this function, you cannot store multiple matrices in 
    different layers of the anndata. 
    
    Input:
    ------
    * rds_file_name [str] --> the name of the inpout file
    
    * path [str] --> the path to the rds file - default: ''
    
    * mtx_name [str] --> matrix/featurespace extracted for the anndata.X - default: 'bmat'
    
    * feature_names [NoneType or str] --> key to extract the variable names (stored as GRanges in
      snap object) if specified - default: 'feature'
    
    * metadata_name [str] --> key to extract the metadata table from the snap object
      (also includes the observation names)- default: 'metaData'
    
    * jaccard [bool] --> extract the jaccard distance if jaccard is True - default: True
    
    * jaccard_key [str] --> If jaccard is True, it will extract the jaccard information
      stored in the snap object under the jaccard_key - default: 'jmat'
    
    * save_all_rds [bool] --> save the entire snap object in the uns section of the Anndata
      object. It takes a lot of memory but can be used to extract the different count matrices
      stored in the snap object - default: False
    
    * verbose [bool] --> print intermediary results and comments if True - default: False
    
    * copy [bool] --> if copy is True, return an Anndata object. - default: True
    
    * save [NoneType or str]--> if save is specified, write the Anndata as .h5ad
       - default: None
    
    
    Output
    ------
    
    * Anndata object --> if copy=True
    
    * Write the anndata object (as h5ad)--> if the output file name is specified as save='output_name'.
        If the path to the output directory is not specified with the output name, the anndata is saved in the current directory.
    
    """
    if path != '':
        file_name = '/'+path.rstrip('/').lstrip('/')+'/'+ rds_file_name.lstrip('/')
    else:
        file_name = rds_file_name.lstrip('/')
    
    tmp_file = extract_snap(rds_file_name=file_name,
                            verbose=verbose)
    
    if copy:
        return(make_Anndata(input_data=tmp_file,
                            mtx_name=mtx_name,
                            feature_names=feature_names,
                            metadata_name=metadata_name, 
                            jaccard=jaccard,
                            jaccard_key=jaccard_key,
                            save_all_rds=save_all_rds,
                            save=save,
                            copy=True))
    else:
        make_Anndata(input_data=tmp_file,
                     mtx_name=mtx_name,
                     feature_names=feature_names,
                     metadata_name=metadata_name, 
                     jaccard=jaccard,
                     jaccard_key=jaccard_key,
                     save_all_rds=save_all_rds,
                     save=save,
                     copy=False)
        