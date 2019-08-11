Installation
------------

Anaconda
~~~~~~~~

If you do not have a working Python 3.5 or 3.6 installation, consider
installing Miniconda (see `Installing Miniconda <http://conda.pydata.org/miniconda.html>`__). Then run::

    conda install seaborn scikit-learn statsmodels numba
    conda install -c conda-forge python-igraph louvain
    conda create -n scanpy python=3.6 scanpy
    
Finally, run::  

    conda install -c annadanese episcanpy



Pull epiScanpy from `PyPI <https://pypi.org/project/episcanpy>`__ (consider
using ``pip3`` to access Python 3)::

    pip install episcanpy
    
    
   
Github
~~~~~~

you can also install epiScanpy directly from Github::

    pip install git+https://github.com/colomemaria/epiScanpy
