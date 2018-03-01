## Overview ##

The goal of this project is to test:

1. When is VPD a driver or reducer of ET?

2. What will drive ET changes with climate? Previous studies (e.g. Jack Scheff's papers) used PET, but the story will be different with ET because of the ecosystem responses. For example, if you include temperature effect on the g<sub>sto</sub> term in Penman-Monteith you will see it cancel a lot of the CC scaling observed in PET. So in reality ET changes as a function of temperature might be more mild, but we need to test this (which is the purpose of this project).

See informal notes in etc/ for more documentation on the genesis of this project.

The general directory structure is:

* src/ - all code
* dat/ - small data files used by code
* doc/ - any papers or formal writing on this project (none yet)
* etc/ - anything that doesn't fit into above, usually notes etc.

### Reproducing research: ###

1) make clean
2) make data
3) make paper

Note that the scripts expect two environmental variables that you will need to set: PLOTS and DATA , which point to arbitrary directories (plots are where misc. plots are saved, and data is where data is loaded from). You will also need to download input data provided by my collaborator Changjie Lin, available at http://www.columbia.edu/~akm2203/data/vpd_data.tar.gz , and extract it into a directory $DATA/changjie/ .  The scripts might also need some directory structure within $PLOTS to successfully save files, so watch for any errors about saving figures; I'll try and check this and update the instructions at some point.

Various latex packages and utilities are required, and below are the details of my python environment used during this project, although note that not all packages listed are required.

Finally, if you run into any issues or bugs, please open an issue.

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