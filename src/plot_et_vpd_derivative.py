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

def concat_dfs(filenames):
  """
  puts all the individual site data into one pdf, and adds a site column to df
  """
  dfs = []
  filenames = glob.glob('%s/changjie/pandas_data/*' % os.environ['DATA'])
  for filename in filenames[:]:
    _df = pd.read_pickle(filename)
    _df['site'] = ''.join(filename.split('/')[-1].split('.')[:-1])
    dfs.append(_df)
  full_df = pd.concat(dfs)
  full_df.to_pickle('%s/changjie/full_pandas.pkl' % os.environ['DATA'])
  return full_df


def make_ax_plot(_ax, var, _df, plot_meta):
  """makes an axis plot"""
  if plot_meta['log'] == 'log':
    var = var/_df['et_obs']
    vmax = 0.8/175.
  elif plot_meta['log'] == 'r_n':
    var = var/(_df['r_n'] - _df['g_flux'])
    vmax = 0.8/175.
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

def scatter_plot(_df, plot_meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 3
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)

  ax1 = fig.add_subplot(nplots, 1, 1)
  var = _df['et_all'] - _df['et']
  plot_meta['cmap'] = 'RdBu'
  plot_meta['delta'] = 'full'
  make_ax_plot(ax1, var, _df, plot_meta)

  ax2 = fig.add_subplot(nplots, 1, 2)
  var = _df['et_leaf'] - _df['et']
  plot_meta['cmap'] = 'RdBu'
  plot_meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, _df, plot_meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  var = _df['et_atm'] - _df['et']
  plot_meta['delta'] = 'atm'
  make_ax_plot(ax3, var, _df, plot_meta)

  plt.tight_layout()
  # plt.savefig('%s/climate_et/site_plots/%s_%s_vpd_debug.png'\
  #             % (os.environ['PLOTS'], str(_df['pft'][0]), plot_meta['site'],))
  fname = '%s/climate_et/%s_plots/%s_%s_%s_%s.png'\
          % (os.environ['PLOTS'], plot_meta['folder_label'],\
             str(_df['pft'][0]), plot_meta['label'],\
             plot_meta['log'], plot_meta['x_axis'])
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


df = pd.read_pickle('%s/changjie/full_pandas.pkl' % os.environ['DATA'])
print(df.shape)
min_et = 0. # W/m2
df = df.loc[(df.et_obs > min_et) & (df.et > min_et), :]
# print(df.shape)
# min_diff = 50. # W/m2
# df = df.loc[(np.absolute(df.et_obs - df.et) < min_diff), :]
# print(df.shape)

plt.close('all')
plot_meta = {}
for x_axis in ['vpd', 'rh']:
  for log in ['log', 'r_n', '']:
    plot_meta['label'] = 'full_ds' 
    plot_meta['folder_label'] = 'full_ds'
    plot_meta['x_axis'] = x_axis
    plot_meta['log'] = log
    # plot_wrapper(df, plot_meta)
    plot_meta['folder_label'] = 'pft'
    df.groupby('pft').apply(plot_wrapper, plot_meta)
    # plot_meta['folder_label'] = 'site'
    # df.groupby('site').apply(plot_wrapper, plot_meta)

# os.system('convert +append %s/climate_et/pft_plots/*false_rh.png '\
#           '%s/climate_et/false_rh_full.png'\
#           % (os.environ['PLOTS'], os.environ['PLOTS']))

# os.system('convert +append %s/climate_et/pft_plots/*true_rh.png '\
#           '%s/climate_et/true_rh_full.png'\
#           % (os.environ['PLOTS'], os.environ['PLOTS']))

# os.system('convert +append %s/climate_et/pft_plots/*true_vpd.png '\
#           '%s/climate_et/true_vpd_full.png'\
#           % (os.environ['PLOTS'], os.environ['PLOTS']))

# os.system('convert +append %s/climate_et/pft_plots/*false_vpd.png '\
#           '%s/climate_et/false_vpd_full.png'\
#           % (os.environ['PLOTS'], os.environ['PLOTS']))
