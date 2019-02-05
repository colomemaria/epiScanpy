Installation
------------

Anaconda
~~~~~~~~

If you do not have a working Python 3.5 or 3.6 installation, consider
installing Miniconda_ (see `Installing Miniconda`_). Then run::

    conda install seaborn scikit-learn statsmodels numba
    conda install -c conda-forge python-igraph louvain
    conda create -n scanpy python=3.6 scanpy



Pull Scanpy from `PyPI <https://pypi.org/project/episcanpy>`__ (consider
using ``pip3`` to access Python 3)::

    pip install episcanpy