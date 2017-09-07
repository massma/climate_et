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
# import matplotlib as mpl
# mpl.rcParams.update(mpl.rcParamsDefault)

resource.setrlimit(resource.RLIMIT_AS, (4000e6, 4000e6))

def test_savefig(fname):
  """tries to save a figure and makes a folder if it doesn't exist"""
  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  return

def split_df(_df):
  """
  splits a pandas df into the atmos, canopy, and data components used
  by various functions on this project. Sometimes it's better or worse
  to have thema ll together.
  """
  atmos = _df.iloc[:, :18].copy()
  canopy = _df.iloc[:, 18:27].copy()
  data = _df.iloc[:, 27:].copy()
  return atmos, canopy, data

def unpack_df(filename):
  """unpacks a dataframe into atmos, canopy, and data components"""
  _df = pd.read_pickle(filename)
  atmos, canopy, data = split_df(_df)
  site = ''.join(filename.split('/')[-1].split('.')[:-1])
  return atmos, canopy, data, site

def concat_dfs(folder='pandas_data_v2', fname='full_pandas_v2'):
  """
  puts all the individual site data into one pdf, and adds a site column to df
  """
  dfs = []
  filenames = glob.glob('%s/changjie/%s/*'\
                        % (os.environ['DATA'], folder))
  for filename in filenames[:]:
    _df = pd.read_pickle(filename)
    _df['site'] = ''.join(filename.split('/')[-1].split('.')[:-1])
    dfs.append(_df)
  full_df = pd.concat(dfs)
  full_df.to_pickle('%s/changjie/%s.pkl'\
                    % (os.environ['DATA'], fname))
  return full_df


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
      vmax = 0.8 # *_df.et_obs.mean()
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
    ax.append(fig.add_subplot(nplots, 1, i+1))
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
  test_savefig(fname)
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
    var = meta['std']*_df['scaling']*(_df['vpd_leaf'])
  else:
    var = _df['scaling']*(_df['vpd_leaf'])
  meta['cmap'] = 'RdBu'
  meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, _df, meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  if meta['var'] == 'd_et_vpd_std':
    var = meta['std']*_df['scaling']*(_df['vpd_atm'])
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

def clean_df(_df, var='lai'):
  """remove unphysical LAI values from a df"""
  out = _df.loc[((_df[var] > 0.1) & (_df[var] < 100.)), :]
  return out

def site_clean(_df, var='lai'):
  """this will remove some percentile of data"""
  out = _df.loc[((_df[var] < _df[var].quantile(q=0.95)) & \
                (_df[var] > _df[var].quantile(q=0.05))), :]
  return out

def histogram(_df, meta):
  """takes a groupby _df and makes histogram plots"""
  fig = plt.figure()
  ax = fig.add_subplot(111)
  sns.distplot(_df[meta['var']], ax=ax)
  ax.set_xlabel(meta['var'])
  if meta['folder_label'] == 'site':
    outname = '%s/%s_%s.png' %\
              (meta['folder'], _df.site.iloc[0], meta['var'])
  elif meta['folder_label'] == 'pft':
    outname = '%s/%s_%s.png' %\
             nn (meta['folder'], _df.pft.iloc[0], meta['var'])
  plt.savefig('%s/climate_et/%s' % (os.environ['PLOTS'], outname))
  return

def test_trend(_df, meta):
  """
  plots the trend of the lai parameter
  to make sure it is independent of vpd
  """
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
  if meta['full_ds']:
    test_savefig('%s/climate_et/scatters/%s_%s.png'\
                % (os.environ['PLOTS'], meta['x_var'], meta['y_var']))
  elif meta['group'] == 'site':
    test_savefig('%s/climate_et/scatters/%s_%s_site/%s.png'\
                 % (os.environ['PLOTS'], meta['x_var'],\
                    meta['y_var'],  _df['site'].iloc[0]))
  else:
    test_savefig('%s/climate_et/scatters/%s_%s/%s.png'\
                % (os.environ['PLOTS'], meta['x_var'],\
                   meta['y_var'],  _df['pft'].iloc[0]))
  return

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
  test_savefig('%s/climate_et/plant_height/%s.png'\
               % (os.environ['PLOTS'], _df.site.iloc[0]))
  return

plt.close('all')

reload_data = False
master_plot = True
if reload_data:
  concat_dfs(folder='pandas_data_lai', fname='full_pandas_lai')
  df = pd.read_pickle('%s/changjie/full_pandas_lai.pkl' % os.environ['DATA'])
  meta = {}
  meta['folder_label'] = 'site'
  meta['folder'] = 'hist_plots'
  meta['var'] = 'lai_gpp'
  print(df.shape)
  df = df.groupby('site').apply(site_clean)
  print(df.shape)
  df = clean_df(df)
  df = clean_df(df, var='lai_gpp')
  # test = df.groupby('site').apply(site_clean, 'lai_gpp')
  # test = clean_df(test, var='lai_gpp')
  df.to_pickle('%s/changjie/full_pandas_lai_clean.pkl' % os.environ['DATA'])
  print(df.shapnne)
  #df.groupby('site').apply(histogram, meta)
  # histogram(df, meta)
  # meta['var'] = 'lai'
  # histogram(df, meta)

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])


# # df.groupby('site').apply(plot_height)
# meta = {}
# meta['xlim'] = None
# meta['ylim'] = None
# meta['plot_type'] = '' #'simple'
# meta['x_var'] = 'swc'
# meta['y_var'] = 'vpd'
# scatter_wrapper(df, meta)
# meta['group'] = 'site'
# df.groupby('site').apply(test_trend, meta)

# meta = {}
# meta['xlim'] = None
# meta['ylim'] = None
# meta['plot_type'] = '' #'simple'
# meta['x_var'] = 'swc'
# meta['y_var'] = 'et_obs'
# scatter_wrapper(df, meta)
# meta['group'] = 'site'
# df.groupby('site').apply(test_trend, meta)


# meta = {}
# meta['x_var'] = 'vpd'
# meta['y_var'] = 'lai'
# meta['xlim'] = (0., 5000.)
# meta['ylim'] = (0.1, 2.)
# for meta['y_var'] in ['lai', 'lai_gpp']:
#   print(meta['y_var'])
#   scatter_wrapper(df, meta)
# meta['xlim'] = None
# meta['ylim'] = None
# meta['x_var'] = 'lai'
# meta['y_var'] = 'lai_gpp'
# scatter_wrapper(df, meta)

if master_plot:
  df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
  df['d_et_atm'] = df['scaling']*df['vpd_atm']

  meta = {}
  meta['vmax'] = None
  meta['var'] = ''
  meta['std'] = df['et_obs'].std()
  for meta['var'], meta['vmax'] in zip(['', 'd_et_vpd_std'], [None, 140.]):
    for x_axis in ['rh']:#'vpd'
      for log in ['scaling', '']:#'log'
        meta['label'] = 'full_ds'
        meta['folder_label'] = 'full_ds'
        meta['x_axis'] = x_axis
        meta['log'] = log
        plot_wrapper(df, meta)
        meta['folder_label'] = 'pft'
        df.groupby('pft').apply(plot_wrapper, meta)
        # meta['folder_label'] = 'site'
        # df.groupby('site').apply(plot_wrapper, meta)

  os.system('convert +append %s/climate_et/pft__rh_plots/*.png '\
            '%s/climate_et/rh.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  os.system('convert +append %s/climate_et/pft_scaling_rh_plots/*.png '\
            '%s/climate_et/rh_scaling.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  os.system('convert +append %s/climate_et/'\
            'd_et_vpd_stdpft_scaling_rh_plots/*.png '\
            '%s/climate_et/d_et_vpd_std_rh_scaling.png'\
            % (os.environ['PLOTS'], os.environ['PLOTS']))

  meta = {}
  meta['log'] = ''
  meta['x_axis'] = 'rh'
  var_lim = {'d_et' : None, 'd_gpp' : None, 'd_wue' : 5.e-5,\
             'd_et_leaf' : None, 'd_et_atm' : None,\
             'vpd_leaf' : 10., 'vpd_atm' : 10., 'scaling' : 0.35}
  var_lim = {'d_et_vpd_std' : 140.,\
             'd_gpp_vpd_std' : 170., 'd_wue_vpd_std' : 0.06}
  for var in var_lim:
    meta['var'] = var
    print(meta['var'])
    meta['vmax'] = var_lim[var]
    meta['folder_label'] = 'full_ds_swc'
    soil_moisture_scatter(df, meta)
    meta['folder_label'] = 'pft_swc'
    df.groupby('pft').apply(soil_moisture_scatter, meta)
    print('working on %s' % var)
    os.system('convert +append %s/climate_et/pft_swc_%s__rh_plots/*.png '\
              '%s/climate_et/swc_%s_rh.png'\
              % (os.environ['PLOTS'], var, os.environ['PLOTS'], var))


