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
# import matplotlib as mpl
# mpl.rcParams.update(mpl.rcParamsDefault)

resource.setrlimit(resource.RLIMIT_AS, (4000e6, 4000e6))

def split_df(_df):
  """
  splits a pandas df into the atmos, canopy, and data components used
  by various functions on this project. Sometimes it's better or worse
  to have thema ll together.
  """
  atmos = _df.iloc[:, :18]# .copy()
  canopy = _df.iloc[:, 18:27]# .copy()
  data = _df.iloc[:, 27:]# .copy()
  return atmos, canopy, data

def unpack_df(filename):
  """unpacks a dataframe into atmos, canopy, and data components"""
  _df = pd.read_pickle(filename)
  atmos, canopy, data = split_df(_df)
  site = ''.join(filename.split('/')[-1].split('.')[:-1])
  return atmos, canopy, data, site



def make_ax_plot(_ax, var, _df, meta):
  """makes an axis plot"""
  if meta['vmax'] is None:
    if meta['log'] == 'log':
      var = var/_df['et_obs']
      vmax = 0.8/175.
    elif meta['log'] == 'scaling':
      var = var/_df['scaling']
      vmax = 5.
    else:
      vmax = 0.7 # *_df.et_obs.mean()
  else:
    vmax = meta['vmax']
  nstd = 1.0
  print(_df['pft'][0])
  # print(var.std())
  print('mean', var.mean())
  vmin = -vmax
  if meta['cmap'] == 'viridis':
    vmax = var.mean() + 0.8 # 5.*var.mean() #nstd*var.std()
    vmin = var.mean() - 0.8 # 5.*var.mean() #nstd*var.std()

  color = _ax.scatter(_df[meta['x_axis']], _df['t_a'], c=var, alpha=0.5,\
                      s=meta['size'], cmap=meta['cmap'],\
                      vmin=vmin, vmax=vmax)
  if (meta['x_axis'] == 'vpd'):
    t_a = np.linspace(_df['t_a'].min(),_df['t_a'].max(), 200.)
    test = met.vapor_pres(t_a)*100.*(1. - 0.90)
    _ax.plot(test, t_a, 'k-')
    test = met.vapor_pres(t_a)*100.*(1. - 0.2)
    _ax.plot(test, t_a, 'k-')
  _ax.set_xlabel(meta['x_axis'])
  _ax.set_ylabel('T')
  _ax.set_title('Analysis: %s,  PFT: %s; %s VPD Changing'\
                % (meta['label'], str(_df['pft'][0]),\
                   meta['delta']))
  cbar = plt.colorbar(color)
  cbar.set_label(r'$\frac{\partial ET}{\partial VPD_{%s}}$'\
                 ' ($W m^{-2}$  $Pa^{-1}$)' % meta['delta'])
  return

def soil_moisture_scatter(_df, meta):
  """
  plots scatter_plots, but going down is different
  percentiels of soil moisture insetad of leaf VPD, etc.
  """
  fig = plt.figure()
  meta['size'] = 16
  percentiles = [ .25, .50, .75, 1.00]
  nplots = len(percentiles)

  fig.set_figheight(fig.get_figheight()*nplots)
  ax = []

  for i, percentile in enumerate(percentiles):
    print('test')
    ax.append(fig.add_subplot(nplots, 1, i+1))
    print(percentile-0.25)
    print(percentile)
    _data = _df.loc[(_df['swc'] > _df.swc.quantile(q=percentile-0.25)) &\
                  (_df['swc'] <= _df.swc.quantile(q=percentile)), :]
    print(_data.shape)
    meta['label'] = 'SWC %d-%d'\
                         % (int((percentile-0.25)*100.), int(percentile*100.))
    meta['cmap'] = 'RdBu'
    meta['delta'] = 'full'
    make_ax_plot(ax[-1], _data[meta['var']], _data, meta)

  plt.tight_layout()
  fname = '%s/climate_et/%s_%s_%s_%s_plots/%s_%s.png'\
          % (os.environ['PLOTS'], meta['folder_label'],\
             meta['var'], meta['log'], meta['x_axis'],
             str(_df['pft'][0]), meta['label'])
  util.test_savefig(fname)
  plt.show(block=False)
  return

def scatter_plot(_df, meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 3
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)
  meta['size'] = 1
  ax1 = fig.add_subplot(nplots, 1, 1)
  # var = _df['et_all'] - _df['et']
  if meta['var'] == 'numeric':
    print('NUMERIC!!!')
    var = _df['et_all']-_df['et']
  elif meta['var'] == 'd_et_vpd_std':
    var = _df['d_et_vpd_std']
  else:
    var = _df['scaling']*(_df['vpd_atm'] + _df['vpd_leaf'])
  meta['log'] = log
  meta['cmap'] = 'RdBu'
  meta['delta'] = 'full'
  make_ax_plot(ax1, var, _df, meta)

  ax2 = fig.add_subplot(nplots, 1, 2)
  if meta['var'] == 'd_et_vpd_std':
    var = _df['d_et_vpd_std_leaf']
  else:
    var = _df['scaling']*(_df['vpd_leaf'])
  meta['cmap'] = 'RdBu'
  meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, _df, meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  if meta['var'] == 'd_et_vpd_std':
    var = _df['d_et_vpd_std_atm']
  else:
    var = _df['scaling']*(_df['vpd_atm'])
  meta['delta'] = 'atm'
  make_ax_plot(ax3, var, _df, meta)

  plt.tight_layout()
  # plt.savefig('%s/climate_et/site_plots/%s_%s_vpd_debug.png'\
  #             % (os.environ['PLOTS'], str(_df['pft'][0]), meta['site'],))
  fname = '%s/climate_et/%s%s_%s_%s_plots/%s_%s.png'\
          % (os.environ['PLOTS'], meta['var'], meta['folder_label'],\
             meta['log'], meta['x_axis'],
             str(_df['pft'][0]), meta['label'])
  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

def plot_wrapper(_df, meta):
  """takes a groupby _df and parses it to plot"""
  print(_df.shape)
  meta['label'] = 'pft'
  if meta['folder_label'] == 'site':
    meta['label'] = str(_df['site'][0])
  scatter_plot(_df, meta)
  return

def gaussian(x, mu, sigma):
  """returns P(x | mu, sigma)"""
  return 1./np.sqrt(2.*np.pi*sigma**2)\
    *np.exp(-(x-mu)**2/(2.*np.sigma**2))

def histogram(_df, meta):
  """takes a groupby _df and makes histogram plots"""
  var = _df[meta['var']]
  fig = plt.figure()
  ax = fig.add_subplot(111)
  sns.distplot(var, ax=ax)
  ax.set_xlabel(meta['var'])
  lims = ax.get_xlim()
  x = np.linspace(lims[0],lims[1], 200)
  ax.plot(x, var.mean(), var.std())
  ax.set_title('pft: %s, site: %s' % (_df.pft.iloc[0], _df.site.iloc[0]))
  if meta['folder_label'] == 'site':
    outname = '%s/%s_%s_%s.png' %\
              (meta['folder'], _df.pft.iloc[0],  _df.site.iloc[0], meta['var'])
  elif meta['folder_label'] == 'pft':
    outname = '%s/%s_%s.png' %\
              (meta['folder'], _df.pft.iloc[0], meta['var'])
  util.test_savefig('%s/climate_et/histogram/%s'\
                    % (os.environ['PLOTS'], outname))
  return

def test_trend(_df, meta, fig=None):
  """
  plots the trend of the lai parameter
  to make sure it is independent of vpd
  """
  if fig is None:
    fig = plt.figure()
  if meta['plot_type'] == 'simple':
    ax = fig.add_subplot(111)
    ax.scatter(_df[meta['x_var']], _df[meta['y_var']], s=8)
    ax.set_xlabel(meta['x_var'])
    ax.set_ylabel(meta['y_var'])
    ax.set_title('spearmanr = %f'\
                 % spearmanr(_df[meta['x_var']],\
                             _df[meta['y_var']]).correlation)
    # ax.set_xlim([0.,1.])
    # ax.set_ylim([0.,1.])
  else:
    g = sns.jointplot(x=_df[meta['x_var']], y=_df[meta['y_var']], kind='hex',\
                      xlim=meta['xlim'], ylim=meta['ylim'], stat_func=spearmanr)
    g.set_axis_labels(meta['x_var'],meta['y_var'])
  plt.title('pft: %s' % _df['pft'].iloc[0])
  if meta['full_ds']:
    util.test_savefig('%s/climate_et/scatters/%s_%s.png'\
                % (os.environ['PLOTS'], meta['x_var'], meta['y_var']))
  elif meta['group'] == 'site':
    util.test_savefig('%s/climate_et/scatters/%s_%s_site/%s.png'\
                 % (os.environ['PLOTS'], meta['x_var'],\
                    meta['y_var'],  _df['site'].iloc[0]))
  else:
    util.test_savefig('%s/climate_et/scatters/%s_%s/%s.png'\
                % (os.environ['PLOTS'], meta['x_var'],\
                   meta['y_var'],  _df['pft'].iloc[0]))
  return spearmanr(_df[meta['x_var']], _df[meta['y_var']]).correlation

def scatter_wrapper(df, meta):
  """just saves line space my wrapping the steps I always take"""
  if 'plot_type' not in meta:
    meta['plot_type'] = 'nan'
  if 'group' not in meta:
    meta['group'] = ''
  meta['full_ds'] = True
  test_trend(df, meta)
  meta['full_ds'] = False
  df.groupby('pft').apply(test_trend, meta)
  return


def plot_height(_df):
  """plots up plant height to make sure it varies"""
  plt.figure()
  plt.plot(np.linspace(0.,100., _df.height.size), _df.height)
  util.test_savefig('%s/climate_et/plant_height/%s.png'\
               % (os.environ['PLOTS'], _df.site.iloc[0]))
  return

def vpd_swc_dependence(_df, meta):
  """at each site plots scatter of swc vs vpd and swc vs et,
  as well as a scatter of det = f(rh, T)"""
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*3)
  axes = [fig.add_subplot(3,1,i+1) for i in range(3)]
  for ax, var in zip(axes, [_df.et_obs, _df.vpd]):
    g = sns.kdeplot(_df.swc, var, shade=True, ax=ax)
    ax.set_title('Spearmanr: %f' % spearmanr(_df.swc, var).correlation)
  meta['log'] = ''
  meta['cmap'] = 'RdBu'
  meta['delta'] = 'full'
  make_ax_plot(axes[2], _df[meta['var']] , _df, meta)
  axes[2].set_title('%s, pft: %s, avg: %f'\
                    % (_df.site.iloc[0],  _df.pft.iloc[0],\
                       _df[meta['var']].mean()))
  plt.tight_layout()
  util.test_savefig('%s/climate_et/scatters/triple_site_%s/%s_%s.png'\
               % (os.environ['PLOTS'], meta['var'],\
                  _df.pft.iloc[0], _df.site.iloc[0]))
  return



plt.close('all')

master_plot = False

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])
df['g_a'] = 1./df['r_a']
df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
df['d_et_atm'] = df['scaling']*df['vpd_atm']

meta = {}
meta['vmax'] = None
meta['var'] = ''
meta['label'] = ''
meta['x_axis'] = 'rh'
meta['log'] = ''
meta['size'] = 8

meta = {}
meta['var'] = 'lai'
meta['folder'] = ''
meta['folder_label'] = site
df.groupby('site').apply(histogram, meta)

# # meta['folder_label'] = 'site'
# for meta['var'], meta['vmax'] in zip(['d_et', 'd_gpp'],\
#                                      [0.3, 0.1]):
#   df.groupby('site').apply(vpd_swc_dependence, meta)

# # df.groupby('site').apply(plot_height)
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
#   scatter_wrapper(df, meta)
#   meta['group'] = 'site'
#   output[name] = df.groupby('site').apply(test_trend, meta)
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
#   scatter_wrapper(df, meta)
#   meta['group'] = 'site'
#   output[name] = df.groupby('site').apply(test_trend, meta)
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
# scatter_wrapper(df, meta)
# meta['group'] = 'site'
# df.groupby('site').apply(test_trend, meta)


# # # meta = {}
# # # meta['x_var'] = 'vpd'
# # # meta['y_var'] = 'lai'
# # # meta['xlim'] = (0., 5000.)
# # # meta['ylim'] = (0.1, 2.)
# # # for meta['y_var'] in ['lai', 'lai_gpp']:
# # #   print(meta['y_var'])
# # #   scatter_wrapper(df, meta)
# # # meta['xlim'] = None
# # # meta['ylim'] = None
# # # meta['x_var'] = 'lai'
# # # meta['y_var'] = 'lai_gpp'
# # # scatter_wrapper(df, meta)

# if master_plot:
#   et_scale = 400.
#   gpp_scale = 400.
#   meta = {}
#   meta['vmax'] = None
#   meta['var'] = ''
#   for meta['var'], meta['vmax'] in zip(['', 'd_et_vpd_std'], [None, et_scale]):
#     for x_axis in ['rh', 'vpd']:
#       for log in ['scaling', '']:#'log'
#         plt.close('all')
#         meta['label'] = 'full_ds'
#         meta['folder_label'] = 'full_ds'
#         meta['x_axis'] = x_axis
#         meta['log'] = log
#         plot_wrapper(df, meta)
#         meta['folder_label'] = 'pft'
#         df.groupby('pft').apply(plot_wrapper, meta)
#         # meta['folder_label'] = 'site'
#         # df.groupby('site').apply(plot_wrapper, meta)

#   os.system('convert +append %s/climate_et/pft__rh_plots/*.png '\
#             '%s/climate_et/rh.png'\
#             % (os.environ['PLOTS'], os.environ['PLOTS']))

#   os.system('convert +append %s/climate_et/pft_scaling_rh_plots/*.png '\
#             '%s/climate_et/rh_scaling.png'\
#             % (os.environ['PLOTS'], os.environ['PLOTS']))

#   os.system('convert +append %s/climate_et/'\
#             'd_et_vpd_stdpft_scaling_rh_plots/*.png '\
#             '%s/climate_et/d_et_vpd_std_rh.png'\
#             % (os.environ['PLOTS'], os.environ['PLOTS']))

#   os.system('convert +append %s/climate_et/'\
#             'd_et_vpd_stdpft_scaling_vpd_plots/*.png '\
#             '%s/climate_et/d_et_vpd_std_vpd.png'\
#             % (os.environ['PLOTS'], os.environ['PLOTS']))

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
#     df.groupby('pft').apply(soil_moisture_scatter, meta)
#     print('working on %s' % var)
#     os.system('convert +append %s/climate_et/pft_swc_%s__rh_plots/*.png '\
#               '%s/climate_et/swc_%s_rh.png'\
#               % (os.environ['PLOTS'], var, os.environ['PLOTS'], var))
