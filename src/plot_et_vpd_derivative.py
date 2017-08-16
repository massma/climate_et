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

def make_ax_plot(_ax, var, atmos, canopy, plot_meta):
  """makes an axis plot"""
  nstd = 2.
  vmax = var.mean() + nstd*var.std()
  vmin = -vmax
  if plot_meta['cmap'] == 'viridis':
    vmax = var.mean() + nstd*var.std()
    vmin = var.mean() - nstd*var.std()
  color = _ax.scatter(atmos['rh'], atmos['t_a'], c=var, alpha=0.2, s=4,\
                         cmap=plot_meta['cmap'], vmin=vmin, vmax=vmax)
  _ax.set_xlabel('RH')
  _ax.set_ylabel('T')
  _ax.set_title('Site: %s; PFT: %s; %s VPD Changing'\
                % (plot_meta['site'], str(canopy['pft'][0]),\
                   plot_meta['delta']))
  cbar = plt.colorbar(color)
  cbar.set_label(r'$\frac{\partial ET}{\partial VPD_{%s}}$'\
                 ' ($W m^{-2}$  $Pa^{-1}$)' % plot_meta['delta'])
  return

def scatter_plot(atmos, canopy, data, site):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 3
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)

  ax1 = fig.add_subplot(nplots, 1, 1)
  var = data['et_all'] - data['et']
  plot_meta = {'cmap' : 'RdBu', 'delta' : 'full', 'site' : site}
  make_ax_plot(ax1, var, atmos, canopy, plot_meta)

  ax2 = fig.add_subplot(nplots, 1, 2)
  var = data['et_leaf'] - data['et']
  plot_meta['cmap'] = 'viridis'
  plot_meta['delta'] = 'leaf'
  make_ax_plot(ax2, var, atmos, canopy, plot_meta)

  ax3 = fig.add_subplot(nplots, 1, 3)
  var = data['et_atm'] - data['et']
  plot_meta['delta'] = 'atm'
  make_ax_plot(ax3, var, atmos, canopy, plot_meta)

  plt.tight_layout()
  plt.savefig('%s/climate_et/site_plots/%s_%s_vpd_debug.png'\
              % (os.environ['PLOTS'], str(canopy['pft'][0]), site,))
  plt.show(block=False)

  return

# filenames = glob.glob('%s/changjie/pandas_data/*' % os.environ['DATA'])
# for filename in filenames[:]:
#   atmos, canopy, data, site = unpack_df(filename)
#   scatter_plot(atmos, canopy, data, site)

