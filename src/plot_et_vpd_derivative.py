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
from scipy.stats import spearmanr
# import matplotlib as mpl
# mpl.rcParams.update(mpl.rcParamsDefault)


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


def make_ax_plot(_ax, var, _df, plot_meta):
  """makes an axis plot"""
  if plot_meta['log'] == 'log':
    var = var/_df['et_obs']
    vmax = 0.8/175.
  elif plot_meta['log'] == 'scaling':
    var = var/_df['scaling']
    vmax = 2000.
  else:
    vmax = 0.8 # *_df.et_obs.mean()
  nstd = 1.0
  print(_df['pft'][0])
  print(var.std())
  print('mean', var.mean())
  # var[var > (var.mean() + nstd*var.std())] = np.nan
  # var[var < (var.mean() - nstd*var.std())] = np.nan
  # vmax = 5.*var.mean() #var.mean() + nstd*var.std()

  vmin = -vmax
  if plot_meta['cmap'] == 'viridis':
    vmax = var.mean() + 0.8 # 5.*var.mean() #nstd*var.std()
    vmin = var.mean() - 0.8 # 5.*var.mean() #nstd*var.std()

  color = _ax.scatter(_df[plot_meta['x_axis']], _df['t_a'], c=var, alpha=0.5,\
                      s=1, cmap=plot_meta['cmap'], vmin=vmin, vmax=vmax)
  if (plot_meta['x_axis'] == 'vpd'):
    t_a = np.linspace(_df['t_a'].min(),_df['t_a'].max(), 200.)
    test = met.vapor_pres(t_a)*100.*(1. - 0.90)
    _ax.plot(test, t_a, 'k-')
    test = met.vapor_pres(t_a)*100.*(1. - 0.2)
    _ax.plot(test, t_a, 'k-')
  _ax.set_xlabel(plot_meta['x_axis'])
  _ax.set_ylabel('T')
  _ax.set_title('Analysis: %s,  PFT: %s; %s VPD Changing'\
                % (plot_meta['label'], str(_df['pft'][0]),\
                   plot_meta['delta']))
  cbar = plt.colorbar(color)
  cbar.set_label(r'$\frac{\partial ET}{\partial VPD_{%s}}$'\
                 ' ($W m^{-2}$  $Pa^{-1}$)' % plot_meta['delta'])
  return

def soil_moisture_scatter(_df, plot_meta):
  """
  plots sam plots as scatter_plot, but going down is different
  percentiels of soil moisture insetad of leaf VPD, etc.
  """
  fig = plt.figure()
  percentiles = [ .25, .50, .75, 1.00]
  nplots = len(percentiles)
  _df['var'] = _df['scaling']*(_df['vpd_atm'] + _df['vpd_leaf'])
  plot_meta['folder_label'] = '%s_%s' % (plot_meta['folder_label'], 'swc')
  fig.set_figheight(fig.get_figheight()*nplots)
  ax = []

  for i, percentile in enumerate(percentiles):
    ax.append(fig.add_subplot(nplots, i+1, 1))
    var = _df['var'].loc[(_df['swc'] > _df.swc.quantile(q=percentile-0.25)) &\
                         (_df['swc'] <= _df.swc.quantile(q=percentile))]
    print(var.shape)
    plot_meta['label'] = 'SWC %d-%d'\
                         % (int((percentile-0.25)*100.), int(percentile*100.))
    plot_meta['log'] = log
    plot_meta['cmap'] = 'RdBu'
    plot_meta['delta'] = 'full'
    make_ax_plot(ax1, var, _df, plot_meta)

  plt.tight_layout()
  fname = '%s/climate_et/%s_%s_%s_plots/%s_%s.png'\
          % (os.environ['PLOTS'], plot_meta['folder_label'],\
             plot_meta['log'], plot_meta['x_axis'],
             str(_df['pft'][0]), plot_meta['label'])
  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

def scatter_plot(_df, plot_meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 3
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)

  ax1 = fig.add_subplot(nplots, 1, 1)
  # var = _df['et_all'] - _df['et']
  if plot_meta['var'] == 'numeric':
    print('NUMERIC!!!')
    var = _df['et_all']-_df['et']
  else:
    var = _df['scaling']*(_df['vpd_atm'] + _df['vpd_leaf'])
  plot_meta['log'] = log
  plot_meta['cmap'] = 'RdBu'
  plot_meta['delta'] = 'full'
  make_ax_plot(ax1, var, _df, plot_meta)

  ax2 = fig.add_subplot(nplots, 1, 2)
  var = _df['scaling']*(_df['vpd_leaf'])
  plot_meta['cmap'] = 'RdBu'
  plot_meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, _df, plot_meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  var = _df['scaling']*(_df['vpd_atm'])
  plot_meta['delta'] = 'atm'
  make_ax_plot(ax3, var, _df, plot_meta)

  plt.tight_layout()
  # plt.savefig('%s/climate_et/site_plots/%s_%s_vpd_debug.png'\
  #             % (os.environ['PLOTS'], str(_df['pft'][0]), plot_meta['site'],))
  fname = '%s/climate_et/%s%s_%s_%s_plots/%s_%s.png'\
          % (os.environ['PLOTS'], plot_meta['var'], plot_meta['folder_label'],\
             plot_meta['log'], plot_meta['x_axis'],
             str(_df['pft'][0]), plot_meta['label'])
  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

def plot_wrapper(_df, plot_meta):
  """takes a groupby _df and parses it to plot"""
  print(_df.shape)
  plot_meta['label'] = 'pft'
  if plot_meta['folder_label'] == 'site':
    plot_meta['label'] = str(_df['site'][0])
  scatter_plot(_df, plot_meta)
  return

def clean_df(_df):
  """remove unphysical LAI values from a df"""
  out = _df.loc[((_df['lai'] > 0.1) & (_df['lai'] < 100.)), :].copy()
  return out

def site_clean(_df):
  """this will remove some percentile of data"""
  out = _df.loc[((_df.lai < _df.lai.quantile(q=0.95)) & \
                (_df.lai > _df.lai.quantile(q=0.05))), :].copy()
  return out

def histogram(_df, plot_meta):
  """takes a groupby _df and makes histogram plots"""
  fig = plt.figure()
  ax = fig.add_subplot(111)
  sns.distplot(_df[plot_meta['var']], ax=ax)
  ax.set_xlabel(plot_meta['var'])
  if plot_meta['folder_label'] == 'site':
    outname = '%s/%s.png' % (plot_meta['folder'], _df.site.iloc[0])
  elif plot_meta['folder_label'] == 'pft':
    outname = '%s/%s.png' % (plot_meta['folder'], _df.pft.iloc[0])
  plt.savefig('%s/climate_et/%s' % (os.environ['PLOTS'], outname))
  return

def test_trend(_df, plot_meta):
  """
  plots the trend of the lai parameter
  to make sure it is independent of vpd
  """
  plt.figure()
  g = sns.jointplot(x=_df['vpd'], y=_df['lai'], kind="hex",\
                    xlim=(0.,5000.), ylim=(0.1,2.), stat_func=spearmanr)
  g.set_axis_labels('vpd','lai')
  if plot_meta['full_ds']:
    plt.savefig('%s/climate_et/lai_vpd.png' % os.environ['PLOTS'])
  else:
    plt.savefig('%s/climate_et/lai_vpd/%s.png'\
                % (os.environ['PLOTS'], _df['pft'].iloc[0]))
  return

# concat_dfs(folder='pandas_data_lai', fname='full_pandas_lai')
# df = pd.read_pickle('%s/changjie/full_pandas_lai.pkl' % os.environ['DATA'])
# plot_meta = {}
# plot_meta['folder_label'] = 'site'
# plot_meta['folder'] = 'hist_plots'
# plot_meta['var'] = 'lai'
# print(df.shape)
# df = df.groupby('site').apply(site_clean)
# print(df.shape)
# df = clean_df(df)
# df.to_pickle('%s/changjie/full_pandas_lai_clean.pkl' % os.environ['DATA'])
# # print(df.shape)
# # df.groupby('site').apply(histogram, plot_meta)



plt.close('all')

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])

test_trend(df, {'full_ds' : True})
df.groupby('pft').apply(test_trend, {'full_ds' : False})
# plot_meta = {}
# plot_meta['var'] = ''
# for x_axis in ['rh', 'vpd']:
#   for log in ['log', 'scaling', '']:
#     plot_meta['label'] = 'full_ds'
#     plot_meta['folder_label'] = 'full_ds'
#     plot_meta['x_axis'] = x_axis
#     plot_meta['log'] = log
#     # plot_meta['var'] = 'numeric'
#     # plot_wrapper(df, plot_meta)
#     plot_meta['folder_label'] = 'pft'
#     df.groupby('pft').apply(plot_wrapper, plot_meta)
#     # plot_meta['folder_label'] = 'site'
#     # df.groupby('site').apply(plot_wrapper, plot_meta)

# # os.system('convert +append %s/climate_et/pft_plots/*false_rh.png '\
# #           '%s/climate_et/false_rh_full.png'\
# #           % (os.environ['PLOTS'], os.environ['PLOTS']))

# # os.system('convert +append %s/climate_et/pft_plots/*true_rh.png '\
# #           '%s/climate_et/true_rh_full.png'\
# #           % (os.environ['PLOTS'], os.environ['PLOTS']))

# # os.system('convert +append %s/climate_et/pft_plots/*true_vpd.png '\
# #           '%s/climate_et/true_vpd_full.png'\
# #           % (os.environ['PLOTS'], os.environ['PLOTS']))

# # os.system('convert +append %s/climate_et/pft_plots/*false_vpd.png '\
# #           '%s/climate_et/false_vpd_full.png'\
# #           % (os.environ['PLOTS'], os.environ['PLOTS']))
