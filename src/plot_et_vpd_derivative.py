#! ~/edward/bin/python
"""
This module plots output from calc_et_vpd_derivative
"""

import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import metcalcs as met
import seaborn as sns
import resource
from scipy.stats import spearmanr
import util
import codebase.plot_tools as plot_tools
# import matplotlib as mpl
# mpl.rcParams.update(mpl.rcParamsDefault)

resource.setrlimit(resource.RLIMIT_AS, (4000e6, 4000e6))

plt.close('all')

master_plot = False

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])
df['g_a'] = 1./df['r_a']
df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
df['d_et_atm'] = df['scaling']*df['vpd_atm']

# meta = {}
# meta['vmax'] = None
# meta['var'] = ''
# meta['label'] = ''
# meta['x_axis'] = 'rh'
# meta['log'] = ''
# meta['size'] = 8

# meta = {}
# meta['var'] = 'lai'
# meta['folder'] = 'site_et'
# meta['folder_label'] = 'site'
# df.groupby('site').apply(plot_tools.histogram, meta)
# meta['folder'] = 'pft_et'
# meta['folder_label'] = 'pft'
# df.groupby('site').apply(plot_tools.histogram, meta)
# meta['folder'] = ''
# meta['folder_label'] = ''
# plot_tools.histogram(df, meta)
# meta['var'] = 'lai_gpp'
# meta['folder'] = 'gpp'
# df.groupby('site').apply(plot_tools.histogram, meta)
# plt.close('all')

# # meta['folder_label'] = 'site'
# for meta['var'], meta['vmax'] in zip(['d_et', 'd_gpp'],\
#                                      [0.3, 0.1]):
#   df.groupby('site').apply(plot_tools.vpd_swc_dependence, meta)

# # df.groupby('site').apply(plot_tools.plot_height)
time = pd.DatetimeIndex(df.time)
df['hour'] = time.hour
df['jd'] = time.dayofyear

# names = ['c_a', 'delta', 'g_a', 'lai', 'vpd', 'jd', 'hour']
# names = ['jd', 'hour']
# output = {}
# for name in names:
#   meta = {}
#   meta['xlim'] = None
#   meta['ylim'] = None
#   meta['plot_type'] = '' #'simple'
#   meta['x_var'] = 'swc'
#   meta['y_var'] = name
#   plot_tools.scatter_wrapper(df, meta)
#   meta['group'] = 'site'
#   output[name] = df.groupby('site').apply(plot_tools.test_trend, meta)
#   plt.close('all')

# for name in output:
#   print('for %s, mean r2 is: %f' % (name, output[name].mean()))
# for name in output:
#   print('for %s, std r2 is: %f' % (name, output[name].std()))

# # names = ['jd', 'hour']
# names = ['c_a', 'delta', 'g_a', 'swc', 'vpd', 'jd', 'hour']
# output = {}
# for name in names:
#   meta = {}
#   meta['xlim'] = None
#   meta['ylim'] = None
#   meta['plot_type'] = '' #'simple'
#   meta['x_var'] = name
#   meta['y_var'] = 'lai'
#   plot_tools.scatter_wrapper(df, meta)
#   meta['group'] = 'site'
#   output[name] = df.groupby('site').apply(plot_tools.test_trend, meta)
#   plt.close('all')

# print('lai')
# for name in output:
#   print('for %s, mean r2 is: %f' % (name, output[name].mean()))
# for name in output:
#   print('for %s, std r2 is: %f' % (name, output[name].std()))

# meta = {}
# meta['xlim'] = None
# meta['ylim'] = None
# meta['plot_type'] = '' #'simple'
# meta['x_var'] = 'swc'
# meta['y_var'] = 'et_obs'
#plot_tools.scatter_wrapper(df, meta)
# meta['group'] = 'site'
# df.groupby('site').apply(plot_tools.test_trend, meta)


# meta = {}
# meta['x_var'] = 'vpd'
# meta['y_var'] = 'lai'
# meta['xlim'] = (0., 5000.)
# meta['ylim'] = (0.1, 2.)
# for meta['y_var'] in ['lai', 'lai_gpp']:
#   print(meta['y_var'])
#  plot_tools.scatter_wrapper(df, meta)
# meta['xlim'] = None
# meta['ylim'] = None
# meta['x_var'] = 'lai'
# meta['y_var'] = 'lai_gpp'
#plot_tools.scatter_wrapper(df, meta)


meta = {}
meta['x_var'] = 'vpd'
meta['y_var'] = 'delta'
meta['xlim'] = None
meta['ylim'] = None
plot_tools.scatter_wrapper(df, meta)

meta = {}
meta['x_var'] = 'rh'
meta['y_var'] = 'delta'
meta['xlim'] = None
meta['ylim'] = None
plot_tools.scatter_wrapper(df, meta)

if master_plot:
  et_scale = 400.
  gpp_scale = 400.
  meta = {}
  meta['vmax'] = None
  meta['var'] = ''
  for meta['var'], meta['vmax'] in zip(['', 'd_et_vpd_std'], [None, et_scale]):
    for x_axis in ['rh', 'vpd']:
      for log in ['scaling', '']:#'log'
        plt.close('all')
        meta['label'] = 'full_ds'
        meta['folder_label'] = 'full_ds'
        meta['x_axis'] = x_axis
        meta['log'] = log
        plot_wrapper(df, meta)
        meta['folder_label'] = 'pft'
        df.groupby('pft').apply(plot_tools.plot_wrapper, meta)
        # meta['folder_label'] = 'site'
        # df.groupby('site').apply(plot_tools.plot_wrapper, meta)

  os.system('convert +append %s/climate_et/pft__rh_plots/*.png '\
            '%s/climate_et/rh.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  os.system('convert +append %s/climate_et/pft_scaling_rh_plots/*.png '\
            '%s/climate_et/rh_scaling.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  os.system('convert +append %s/climate_et/'\
            'd_et_vpd_stdpft_scaling_rh_plots/*.png '\
            '%s/climate_et/d_et_vpd_std_rh.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  os.system('convert +append %s/climate_et/'\
            'd_et_vpd_stdpft_scaling_vpd_plots/*.png '\
            '%s/climate_et/d_et_vpd_std_vpd.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

#   meta = {}
#   meta['log'] = ''
#   meta['x_axis'] = 'rh'
#   var_lim = {'d_et' : None, 'd_gpp' : None, 'd_wue' : 5.e-5,\
#              'd_et_leaf' : None, 'd_et_atm' : None,\
#              'vpd_leaf' : 10., 'vpd_atm' : 10., 'scaling' : 0.35,\
#              'd_et_vpd_std' : et_scale,\
#              'd_gpp_vpd_std' : gpp_scale, 'd_wue_vpd_std' : 0.06}
#   for var in var_lim:
#     plt.close('all')
#     meta['var'] = var
#     print(meta['var'])
#     meta['vmax'] = var_lim[var]
#     meta['folder_label'] = 'full_ds_swc'
#     soil_moisture_scatter(df, meta)
#     meta['folder_label'] = 'pft_swc'
#     df.groupby('pft').apply(plot_tools.soil_moisture_scatter, meta)
#     print('working on %s' % var)
#     os.system('convert +append %s/climate_et/pft_swc_%s__rh_plots/*.png '\
#               '%s/climate_et/swc_%s_rh.png'\
#               % (os.environ['PLOTS'], var, os.environ['PLOTS'], var))

