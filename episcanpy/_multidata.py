import anndata as ad
import numpy as np
import pandas as pd

import os.path
from os import path
import os
from os import listdir
from os.path import isfile, join

import shutil
import csv
import warnings

class MultiData():
     #"""
     #Object to contain one or more annotated data matrices.
     #"""
    def __init__(self,
                 anndata=None,
                 omic_key=None,
                 barcodes=None,
                 paired=True):
        """
        anndata can be either a unique object or a list of multiple objects. 
        omic_key correspond to the attribute name or the key for the dictionaty
        paired --> paired is still either a bool or a list. if it is a bool, the same thing will be
        remembered for everything
        
        barcodes needs to be a pandas dataframe where the columns correspond to omic_key and the rows
        correspond to the cells
        
        # there will be a problem if the omic key contains a _ in its name
        
        """
        
        ## saev the omic labels
        if isinstance(omic_key, list) is False:
            omic_key = [omic_key]
        self.omic_label = omic_key
           
        ## add the anndata objects
        self.omic = dict()
        if isinstance(anndata, list) is False:
            anndata = [anndata]
        index = 0
        for key in omic_key:
            self.omic[key] = anndata[index]
            index += 1
        
        ## add the paired objects
        if isinstance(paired, bool):
            self.paired = [paired]*len(omic_key)
        else:
            self.paired = paired
        
        self.barcodes = []
        self.uns = {}
    
        
    def __repr__(self) -> str:
        """
        """
        descr = f"MultiData object contains {len(self.omic)} Anndata with the follogwing keys {self.omic.keys()}"
        #descr = f"AnnData object with n_obs × n_vars = {n_obs} × {n_vars}{backed_at}"
        for attr in ["omic", "paired", "omic_label", "barcodes", "uns"]:
            keys = getattr(self, attr)
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return(descr)
        
    
    #### write a multidata folder function
    
    def _write(self, folder_name):
        """
        """
        for attr in ['omic', 'paired', 'omic_label']:
            
            data = getattr(multi, attr) # extract the data in the attribute
            output_file_name = folder_name+"/"+folder_name.split('/')[-1]+"_"+attr
            
            if attr == 'omic': # save the multiple anndata
                for omic_layer in data.keys():
                    data[omic_layer].write(output_file_name+"_"+omic_layer+".h5ad")
                    
                # save in which order do we need to read the different anndata   
                with open(output_file_name+"_omic_keys"+'.csv','w') as f:
                    writer = csv.writer(f)
                    writer.writerow(data.keys())
                    
            else: # to save the other kind of information contained in the multidata class
                with open(output_file_name+'.csv','w') as f:
                    writer = csv.writer(f)
                    writer.writerow(data)
        
        
    def write(self, folder_name='multidata_object', overwrite=True):
        """
        Create a folder. If the folder already exist, return a warning and add/update h5ad files without removing the old ones
        If you want to overwrite teh folder, there is an overwrite parameter.
        
        folder_name can take the path (absolute or relative) as well as the new directory to create. 
        if a path is not specified, the function creates the new folder in the working directory. 
        """
        # check if path exist:
        if path.exists(folder_name):
            if overwrite==True:
                # check if there is permission to write
                try:
                    shutil.rmtree(folder_name) # delete directory and its content
                    os.mkdir(folder_name) # create the directory again 
                except OSError as e:
                    print("Error: %s : %s" % (folder_name, e.strerror))
            else:
                warnings.warn("""The overwrite paramter is set to false.
                Old element of the multidata will not be updated.\n
                It might create some issues when reading multidata""")
        else:
            os.mkdir(folder_name) # create the output directory
        
        self._write(folder_name)## write function
        
        
    ### add an omic function
    ### matching barcode ? 
    ### merging 2 paired anndata function 

def read_multidata(input_folder):
    """
    read all the files in input_folder and return a multidata_object
    """
    path = input_folder.rstrip('/')+'/'
    # find what is in the folder: 
    filestoload = [f for f in listdir(input_folder) if isfile(join(input_folder, f))]
    #print(filestoload)
    # sort out the different files
    h5ad_files = []
    csv_files = []
    for file_name in filestoload:
        file_name=path+file_name
        #print(file, path)
        if 'omic_keys.csv' in file_name:
            # read the file specifying in which order the omics are saved
            order_file=pd.read_csv(file_name, header=None).iloc[0,:].tolist()
        elif '.h5ad' in file_name:
            h5ad_files.append(file_name)
        else:
            csv_files.append(file_name)
            
    # I need to read the files in order.
    
    # read .h5ad --> save in the omic attribute
    dict_of_anndata = {}
    for adata_key in h5ad_files:
        # there will be a problem if the omic key contains a _ in its name
        dict_of_anndata[adata_key.rstrip('.h5ad').split('_')[-1]] = ad.read(adata_key)
        
    # make the multidata object
    multi = MultiData(anndata=[dict_of_anndata[omic] for omic in order_file], omic_key=order_file)
        
    # read .csv files --> save as attributes
    for file_name in csv_files:
        if "_paired.csv" in file_name:
            multi.paired = pd.read_csv(file_name, header=None).iloc[0,:].tolist()
        elif "omic_label.csv" in file_name:
            multi.paired = pd.read_csv(file_name, header=None).iloc[0,:].tolist()
        else:
            multi.uns[file_name] = pd.read_csv(file_name, header=None).iloc[0,:].tolist()
    
    return(multi)