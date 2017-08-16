#! ~/edward/bin/python
"""
This module plots output from calc_et_vpd_derivative
"""

import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

def make_ax_plot(_ax, var, _df, plot_meta):
  """makes an axis plot"""
  nstd = 2.
  vmax = var.mean() + nstd*var.std()
  vmin = -vmax
  if plot_meta['cmap'] == 'viridis':
    vmax = var.mean() + nstd*var.std()
    vmin = var.mean() - nstd*var.std()
  print(var.shape)
  print(_df['rh'])
  print(_df['t_a'])
  print(plot_meta['cmap'])
  print(vmin)
  print(vmax)
  print(var)
  color = _ax.scatter(_df['rh'], _df['t_a'], c=var, alpha=0.2, s=4,\
                         cmap=plot_meta['cmap'], vmin=vmin, vmax=vmax)
  _ax.set_xlabel('RH')
  _ax.set_ylabel('T')
  _ax.set_title('Site: %s; PFT: %s; %s VPD Changing'\
                % (plot_meta['site'], str(_df['pft'][0]),\
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
  make_ax_plot(ax1, var, df, plot_meta)

  ax2 = fig.add_subplot(nplots, 1, 2)
  var = _df['et_leaf'] - _df['et']
  plot_meta['cmap'] = 'viridis'
  plot_meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, df, plot_meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  var = _df['et_atm'] - _df['et']
  plot_meta['delta'] = 'atm'
  make_ax_plot(ax3, var, df, plot_meta)

  plt.tight_layout()
  # plt.savefig('%s/climate_et/site_plots/%s_%s_vpd_debug.png'\
  #             % (os.environ['PLOTS'], str(_df['pft'][0]), plot_meta['site'],))
  plt.savefig('%s/climate_et/%s_plots/%s_vpd_debug.png'\
              % (os.environ['PLOTS'], plot_meta['label'], str(_df['pft'][0])))

  plt.show(block=False)
  return

def plot_wrapper(_df, *args):
  """takes a groupby _df and parses it to plot"""
  plot_meta = {}
  plot_meta['label'] = args
  scatter_plot(_df, plot_meta)
  return

# filenames = glob.glob('%s/changjie/pandas_data/*' % os.environ['DATA'])
# for filename in filenames[:]:
#   atmos, canopy, data, site = unpack_df(filename)
#   scatter_plot(atmos, canopy, data, site)

df = pd.read_pickle('%s/changjie/full_pandas.pkl' % os.environ['DATA'])
df.groupby('pft').apply(plot_wrapper, 'pft')
