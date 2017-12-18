#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *

### FIGURE 2 ###

meta = {}
meta['var'] = 'lai'
meta['folder'] = ''
meta['folder_label'] = ''
plot_tools.histogram(df, meta)
os.system('cp %s/climate_et/histogram/full_lai.png '\
          '../../doc/paper/fig02.png' % (os.environ['PLOTS']))

### FIGURE 3 ###
# joint distribution of lai and swc
meta = {}
meta['ylim'] = (0., 2.)
meta['xlim'] = (0., 4000.)
meta['plot_type'] = '' #'simple'
meta['x_var'] = 'vpd'
meta['y_var'] = 'lai'
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/vpd_lai.png ../../doc/paper/fig03.png'\
          % (os.environ['PLOTS']))
