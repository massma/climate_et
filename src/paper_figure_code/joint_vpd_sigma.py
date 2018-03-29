#! ~/edward/bin/python
"""
This script makes fig 3 for paper
"""
from shared_functions import *

plt.close('all')

importlib.reload(plot_tools)
### FIGURE 3 ###
# joint distribution of lai and swc
df['sigma'] = df['sigma'].copy()
meta = {}
meta['ylim'] = (0., 2.)
meta['xlim'] = (0., 4000.)
meta['plot_type'] = '' #'simple'
meta['x_var'] = 'vpd'
meta['y_var'] = 'sigma'
meta['x_label'] = 'VPD (Pa)'
meta['y_label'] = '$\sigma$'
meta['pft_plot'] = False
meta['fontsize'] = single_ax_fontsize+3
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/vpd_sigma.pdf '\
          '../../doc/paper/joint_vpd_sigma.pdf'\
          % (os.environ['PLOTS']))

os.system('cp %s/climate_et/scatters/vpd_sigma.pdf '\
          '../../doc/shared_figs/vpd_sigma.pdf'\
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
