Tutorials
=========


------------

Single cell ATAC-seq
--------------------

To get started, we recommend epiScanpy's analysis pipeline for scATAC-seq data from Buenrostro et al. [Buenrostro18]_. , the dataset consist of ~3000cells of human PBMCs. This tutorial focuses on preprocessing, clustering, identification of cell types via known marker genes and trajectory inference. The tutorial can be found `here <https://nbviewer.jupyter.org/github/colomemaria/epiScanpy/blob/master/docs/tutorials/Buenrostro_PBMC_data_processing.html>`__. 



.. image:: pictures/Buenrostro_umap.png
   :width: 200px


If you want to see how to build count matrices from ATAC-seq bam files for different set of annotations (like enhancers).
The tutorial can be found `here <https://nbviewer.jupyter.org/github/colomemaria/epiScanpy/blob/master/docs/tutorials/ATAC_bld_ct_mtx_tutorial.html>`__. 

Soon available, there will be a tutorial providing a function to very quickly build custom count matrices using standard 10x single cell ATAC output. 

.. image:: pictures/umap.png
   :width: 350px
.. image:: pictures/umap_ATACseq_Astrocyte_marker.png
   :width: 250px
.. image:: pictures/heatmap.png
   :width: 600px

An additional tutorial on processing and clustering count matrices from the Cusanovich mouse scATAC-seq atlas [Cusanovich18]_.. `Here <https://nbviewer.jupyter.org/github/colomemaria/epiScanpy/blob/master/docs/tutorials/Cusanovich2018_BoneMarrow_data_processing_diffmap.html>`__.

------------

Single cell DNA methylation
---------------------------

Here you can find a tutorial for the preprocessing, clustering and identification of cell types for single-cell DNA methylation data using the publicly available data from Luo et al. [Luo17]_. 

The first tutorial shows how to build the count matrices for the different feature spaces (windows, promoters) in different cytosine contexts. Here is the  `tutorial  <https://nbviewer.jupyter.org/github/colomemaria/epiScanpy/blob/master/docs/tutorials/bld_count_matrix_methylation_tutorial.html>`__.

To see how process single cell DNA methylation data suing epiScanpy. Here is a `tutorial <https://nbviewer.jupyter.org/github/colomemaria/epiScanpy/blob/master/docs/tutorials/processing_enhancers_mCG.html>`__ processing The data Luo et al. [Luo17]_ mouse frontal cortex dna methylation data. 

.. image:: pictures/umap_markers_hodology_ecker.png
   :width: 600px
.. image:: pictures/umapexcitatory_neurons_promoters.png
   :width: 300px 
.. image:: pictures/umapSatb2_CLUSTER_NORM.png
   :width: 250px  
.. image:: pictures/umapenhancer_CG_Luoetal.png
   :width: 300px
