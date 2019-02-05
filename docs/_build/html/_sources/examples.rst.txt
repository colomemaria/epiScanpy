Tutorials
=========


------------

single cell ATAC
----------------

Currently, for getting started, we recommend `Scanpy's reimplementation <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb>`__ of Seurat's [Satija15]_ clustering tutorial 2700 PBMCs from 10x Genomics, containing preprocessing, clustering and the identification of cell types via known marker genes.

You might compare this to `clustering 68K PBMCs <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170503_zheng17/zheng17.ipynb>`__ as it would be done in Cell Ranger [Zheng17]_.

.. raw:: html

   <img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/filter_genes_dispersion.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/louvain.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/NKG7.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/violin.png" style="width: 100px"><img src="http://falexwolf.de/img/scanpy_usage/170505_seurat/cell_types.png" style="width: 200px">


------------

Single cell methylation
-----------------------

For trajectory inference on complex datasets, we offer several examples `here <https://github.com/theislab/paga>`__. Get started `here <https://nbviewer.jupyter.org/github/theislab/paga/blob/master/blood/paul15/paul15.ipynb>`__ for the following result on hematopoiesis.

.. raw:: html

   <img src="http://www.falexwolf.de/img/paga_paul15.png" style="width: 450px">

You can extend this to multi-resolution analyses of whole animals, such as `here <https://nbviewer.jupyter.org/github/theislab/paga/blob/master/planaria/planaria.ipynb>`__.

.. raw:: html

   <img src="http://www.falexwolf.de/img/paga_planaria.png" style="width: 350px">

The PAGA method behind this is described `here <https://rawgit.com/falexwolf/paga_paper/master/paga.pdf>`__ and can be cited using this `doi <https://doi.org/10.1101/208819>`__. As a reference for simple pseudotime analyses, we provide the diffusion pseudotime analyses of [Haghverdi16]_ for two hematopoiesis datasets: `here <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170502_paul15/paul15.ipynb>`__ for [Paul15]_ and `here <https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170501_moignard15/moignard15.ipynb>`__ for [Moignard15]_.