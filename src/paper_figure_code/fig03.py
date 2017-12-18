#! ~/edward/bin/python
"""
This script makes fig 3 for paper
"""
from shared_functions import *

plt.close('all')

### FIGURE 3 ###
# joint distribution of lai and swc
df['sigma'] = df['sigma'].copy()
meta = {}
meta['ylim'] = (0., 2.)
meta['xlim'] = (0., 4000.)
meta['plot_type'] = '' #'simple'
meta['x_var'] = 'vpd'
meta['y_var'] = 'sigma'
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/vpd_sigma.png '\
          '../../doc/paper/fig03.png'\
          % (os.environ['PLOTS']))

os.system('cp %s/climate_et/scatters/vpd_sigma.png '\
          '../../doc/shared_figs/vpd_sigma.png'\
          % (os.environ['PLOTS']))

# below added for agu, not currently sued in paper but could be
# useful for supplemental material
meta['x_var'] = 'e_s'
meta['y_var'] = 'rh'
meta['ylim'] = None
meta['xlim'] = None
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/e_s_rh.png '\
          '../../doc/shared_figs/e_s_rh.png'\
          % (os.environ['PLOTS']))
