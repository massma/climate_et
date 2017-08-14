#! ~/edward/bin/python
"""
This module plots output from calc_et_vpd_derivative
"""

import time
import glob
import os
import pandas as pd

filenames = glob.glob('%s/changjie/pandas_data/*' % os.environ['DATA'])

def unpack_df(filename):
  """unpacks a dataframe into atmos, canopy, and data components"""
  _df = pd.read_pickle(filename)
  atmos = _df.iloc[:, :18].copy()
  canopy = _df.iloc[:, 18:27].copy()
  data = _df.iloc[:, 27:].copy()
  return atmos, canopy, data
  

def scatter_plot(data):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 3
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)

  ax1 = fig.add_subplot(nplots, 1, 1)
  var = data['et_all'] - data['et']
  #now do scatter plot, maybe using below
  _ax.append(fig.add_subplot(nplots, 1, i+1))
    vmax = result[key].mean() + 3.*result[key].std()
    vmin = -vmax
    if i > 0:
      cmap = 'viridis'
      vmax = None
      vmin = None
    else:
      cmap = 'RdBu'
    color = _ax[i].pcolormesh(atmos['rh'], atmos['t_a'], result[key],\
                             cmap=cmap, vmin=vmin, vmax=vmax)
    _ax[i].set_xlabel('RH')
    _ax[i].set_ylabel('T')
    _ax[i].set_title('PFT: %s; Model: %s; %s VPD Changing'\
                     % (canopy['pft'], canopy['stomatal_model'], str(key)))
    cbar = plt.colorbar(color)
    cbar.set_label(r'$\frac{\partial ET}{\partial VPD_{%s}}$'\
                   ' ($W m^{-2}$  $Pa^{-1}$)' % key)
  plt.tight_layout()
  plt.savefig('%s/climate_et/%s_%s_vpd_debug.png'\
              % (os.environ['PLOTS'], canopy['pft'],\
                 canopy['stomatal_model']))
  plt.show(block=False)

  return
for filename in filenames[:1]:
  atmos, canopy, data = unpack_df(filename)
  
  

