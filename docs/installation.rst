Installation
------------

Anaconda
~~~~~~~~

If you do not have a working Python 3.5 or 3.6 installation, consider
installing Miniconda_ (see `Installing Miniconda`_). Then run::

    conda install seaborn scikit-learn statsmodels numba
    conda install -c conda-forge python-igraph louvain
    conda create -n scanpy python=3.6 scanpy
    
Finally, run::  
    conda install episcanpy



Pull Scanpy from `PyPI <https://pypi.org/project/episcanpy>`__ (consider
using ``pip3`` to access Python 3)::

    pip install episcanpy
    
    
   
Github
~~~~~~

you can also (right know it is the only way available) install
epiScanpy like the following::

    pip install git+https://github.com/colomemaria/epiScanpy



System Requirements
-------------------

Hardware requirements
~~~~~~~~~~~~~~~~~~~~~

``epiScanpy`` package requires only a standard computer with enough RAM to support the in-memory operations.

Software requirements
~~~~~~~~~~~~~~~~~~~~~

### OS Requirements
This package is supported for *macOS* and *Linux*. The package has been tested on the following systems:
+ macOS: Mojave (10.14.1)

Python Dependencies
~~~~~~~~~~~~~~~~~~~
``epiScanpy`` depends on the Python scientific stack.::

  numpy
  pandas
  scanpy
  scipy
  scikit-learn
  seaborn

