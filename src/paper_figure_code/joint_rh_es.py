#! ~/edward/bin/python
"""
This script makes fig 3 for paper
"""
from shared_functions import *

plt.close('all')

importlib.reload(plot_tools)
### FIGURE 3 ###
# joint distribution of lai and swc
df['rh'] = df['rh'].copy()
df['e_s'] = df['e_s'].copy()
meta = {}
meta['xlim'] = (0., 100.0)
meta['ylim'] = (0.0, 7000.0)
meta['plot_type'] = '' #'simple'
meta['x_var'] = 'rh'
meta['y_var'] = 'e_s'
meta['x_label'] = 'Relative Humidity (Percent)'
meta['y_label'] = 'Saturation Vapor Pressure (Pa)'
meta['pft_plot'] = False
meta['fontsize'] = single_ax_fontsize
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/rh_e_s.pdf '\
          '../../doc/paper/supp-figs/0joint_rh_es.pdf'\
          % (os.environ['PLOTS']))

# below added for agu, not currently sued in paper but could be
# useful for supplemental material
# meta={}
# meta['plot_type'] = ''
# meta['x_var'] = 'e_s'
# meta['y_var'] = 'rh'
# meta['ylim'] = None
# meta['xlim'] = None
# plot_tools.scatter_wrapper(df, meta)
# os.system('cp %s/climate_et/scatters/e_s_rh.png '\
#           '../../doc/shared_figs/e_s_rh.png'\
#           % (os.environ['PLOTS']))

# meta={}
# meta['plot_type'] = ''
# meta['x_var'] = 'r_a_corrected'
# meta['y_var'] = 'r_a_uncorrected'
# meta['ylim'] = (0., 100.)
# meta['xlim'] = (0., 100.)
# meta['pft_plot'] = True

# plot_tools.scatter_wrapper(df, meta)
# os.system('cp %s/climate_et/scatters/e_s_rh.png '\
#           '../../doc/shared_figs/e_s_rh.png'\
#           % (os.environ['PLOTS']))
