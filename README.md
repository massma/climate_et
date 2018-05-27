## When does vapor pressure deficit drive or reduce ET? 

This is the repository for reproducing the (currently pre-printed) manuscript:

Massmann, A., P. Gentine and C. Lin. _When does vapor pressure deficit drive or reduce evapotranspiration?_ [arXiv:1805.05444](https://arxiv.org/abs/1805.05444).

The general directory structure is:

* src/ - all code
* dat/ - small data files used by code
* doc/ - any formal papers and presentations for the project 
* etc/ - anything that doesn't fit into above, usually notes etc.

Note that you also need to initialize submodules:

* ```git submodule init && git submodule update```
* ```cd src/codebase/fluxnet_pycite && git submodule init && git submodule update```

### Reproducing research:

1) make clean
2) make data
3) make paper
4) make arxiv (if you want tarball submitted to arxiv)

Note that the scripts expect two environmental variables that you will need to set: PLOTS and DATA , which point to arbitrary directories (plots are where misc. plots are saved, and data is where data is loaded from). The scripts might also need some directory structure within $PLOTS to successfully save files, so watch for any errors about saving figures; I'll try and check this and update the instructions at some point. You will also need to download processed [FLUXNET2015](https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/) data provided by my collaborator Changjie Lin, available at http://www.columbia.edu/~akm2203/data/vpd_data.tar.gz , and extract it into a directory $DATA/changjie/ . It is extremely important that these data are only used for the purposes of reproducing this research, for two good reasons. First, FLUXNET data for new research should be optained from the FLUXNET website to insure that the data policy is followed and FLUXNET is aware of the downloads. Second, only a portion of the data in the tarball are actually used in the analysis, and any of the other variables may not be reliable, quality controlled, or correct.  

Finally, if you run into any issues or bugs, please open an issue or shoot me an email. 

python --version : Python 3.5.3

pip freeze :
astroid==1.5.3
backports.weakref==1.0rc1
basemap==1.0.7
beautifulsoup4==4.5.3
bleach==1.5.0
cairocffi==0.7.2
cffi==1.9.1
chardet==2.3.0
cycler==0.10.0
decorator==4.0.11
edward==1.3.4
epc==0.0.5
flake8==3.2.1
GDAL==2.1.2
h5py==2.7.0
html5lib==0.9999999
httplib2==0.9.2
isort==4.2.15
jedi==0.10.0
joblib==0.10.4.dev0
Keras==2.0.8
lazy-object-proxy==1.3.1
lxml==3.7.1
Markdown==2.2.0
matplotlib==2.0.0
mccabe==0.5.3
mpmath==0.19
netCDF4==1.2.7
nose==1.3.7
numexpr==2.6.1
numpy==1.13.3
observations==0.1.3
pandas==0.21.0
pep8==1.7.0
Pillow==4.0.0
pkg-resources==0.0.0
plotly==1.13.0
ply==3.9
protobuf==3.3.0
py==1.4.32
pycodestyle==2.2.0
pycparser==2.17
pycurl==7.43.0
pyflakes==1.3.0
pygobject==3.22.0
pylint==1.7.2
pyparsing==2.1.10
pytest==3.0.6
python-apt==1.4.0b3
python-dateutil==2.6.1
python-debian==0.1.30
python-debianbts==2.6.1
pytz==2017.3
PyYAML==3.12
reportbug==7.1.7
requests==2.12.4
scikit-learn==0.18
scipy==0.18.1
seaborn==0.7.1
sexpdata==0.0.3
simplejson==3.10.0
six==1.11.0
sympy==1.0
tables==3.3.0
tensorflow==1.2.0
urllib3==1.19.1
virtualenv==15.1.0
webencodings==0.5
Werkzeug==0.12.2
wrapt==1.10.10
xarray==0.9.6