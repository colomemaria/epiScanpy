import sys, os
if sys.version_info < (3,):
    sys.exit('epiScanpy requires Python >= 3.5')
from pathlib import Path
from setuptools import setup, find_packages
import versioneer

try:
    from episcanpy import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ''

# extract version
#with open(os.path.join(os.path.dirname(__file__), "episcanpy", "version.py")) as f:
#    version = f.read().split("\n")[0].split("=")[-1].strip(' ').strip('"')

setup(
    name='episcanpy',
    version=versioneer.get_version(),
    #version=version,
    cmdclass=versioneer.get_cmdclass(),
    description='Epigenomics Single-Cell Analysis in Python.',
    long_description=Path('README.rst').read_text('utf-8'),
    #long_description=read_file('README.rst'),
    url='http://github.com/colomemaria/epiScanpy',
    #download_url='https://github.com/colomemaria/epiScanpy/archive/v0.1-beta.tar.gz',
    author=__author__,
    author_email=__email__,
    license='BSD',
    python_requires='>=3.5',
    install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    extras_require=dict(
        louvain=['python-igraph', 'louvain>=0.6'],
        leiden=['python-igraph', 'leidenalg'],
        bbknn=['bbknn'],
        combat=['patsy'],
        doc=['sphinx', 'sphinx_rtd_theme', 'sphinx_autodoc_typehints', 'scanpydoc'],
        test=['pytest'],
    ),
    packages=find_packages(),
    # `package_data` does NOT work for source distributions!!!
    # you also need MANIFTEST.in
    # https://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
