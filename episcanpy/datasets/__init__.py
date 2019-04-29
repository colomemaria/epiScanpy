import anndata as ad
import scanpy.api as sc

def Luo17_raw():
    """Methylation summaries for 5 cells from the mouse prefrontal cortex [Luo17]_.
    
    Returns
    -------
    Methylation summaries name files.
    """
    [...]
    return(methylation_summaries)

def Luo17promoters_raw():
    """CG Methylation levels at promoters. Mouse prefrontal cortex neurons [Luo17]_.
    
    Returns
    -------
    Annotated data matrix (not filtered).
    """ 
    ## Here I put the raw count matrix. Not the annotated/filtered one 
  [...]
  return(count_matrix)
  
def Luo17promoters():
    """CG Methylation levels at promoters. Mouse prefrontal cortex neurons [Luo17]_.
    
    Returns
    -------
    Annotated data matrix.
    """ 
 
 # I will need to add the hematopoitic count matrix from Cusanovich and the peak count matrix for brain data
 # Also Cusanovich et al. 2018
