#! ~/edward/bin/python
"""
This script makes fig 3 for paper
"""
from shared_functions import *

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
