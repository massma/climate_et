#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *

### FIGURE 2 ###

meta = {}
meta['var'] = 'sigma'
meta['folder'] = ''
meta['folder_label'] = ''
plot_tools.histogram(df, meta)
os.system('cp %s/climate_et/histogram/full_sigma.png '\
          '../../doc/paper/fig02.png' % (os.environ['PLOTS']))

