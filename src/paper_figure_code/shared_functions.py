#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
import os
import importlib
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import codebase.plot_tools as plot_tools
import util
import codebase.data_calc as d_calc
import codebase.data_prep as d_prep
import pickle

mpl.use('Pdf')
mpl.rcParams.update(mpl.rcParamsDefault)
small_ax_params = mpl.rcParamsDefault
# increase tick label size
label_size = 14
small_ax_params['xtick.labelsize'] = label_size
small_ax_params['ytick.labelsize'] = label_size


importlib.reload(util)
importlib.reload(plot_tools)
importlib.reload(d_calc)

# grab sites actually used in analysis
with open('%s/changjie/diagnosed_data.pkl' % os.environ['DATA'],\
          mode='rb') as file:
  dfs = pickle.load(file)
df = dfs['full']
median_df = dfs['median']
mean_df = dfs['mean']

def et_min_vpd(_df, uwue=None):
  """calculates theoretical vpd_crit as functoin of -df and lai"""
  if uwue is None:
    uwue = _df.uwue
  c3 = d_calc.CP/_df.r_moist
  c1 = _df.gamma*_df.c_a/(d_calc.R_STAR*1.6*uwue)
  c2 = _df.g1
  sqrt_vpd = (c1 + np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3)
  try:
    sqrt_vpd[sqrt_vpd < 0.] = np.nan
  except TypeError:
    if sqrt_vpd < 0.:
      sqrt_vpd = np.nan
  return sqrt_vpd**2

def get_pft(_df):
  return _df['pft'].iloc[0]

def d_et(_df):
  """returns d et/d vpd given _df"""
  return d_calc.scaling(_df)*d_calc.sign(_df)

name_dict = {'CRO': 'Crops (CRO)',\
             'DBF': 'Deciduous Forest (DBF)',
             'EBF': 'Evergreen Broadleaf Forest (EBF)',
             'ENF': 'Evergreen Needleleaf Forest (ENF)',
             'GRA': 'Grass (GRA)',
             'CSH': 'Closed Shrub (CSH)',
             'WSA': 'Woody Savannah (WSA)',
             'MF': 'Mixed Forest (MF)',
             'SAV': 'Savannah (SAV)'}

fontsize=16

def custom_ylabel(_df, ax, label):
  """labels y label for axis ordering"""
  if (_df.pft.iloc[0] == 'ENF') | (_df.pft.iloc[0] == 'DBF')\
     | (_df.pft.iloc[0] == 'CRO'):
    ax.set_ylabel(label, fontsize=fontsize)
  return

def custom_xlabel(_df, ax, label):
  """labels x label for axis ordering"""
  if (_df.pft.iloc[0] == 'CRO') | (_df.pft.iloc[0] == 'GRA')\
     | (_df.pft.iloc[0] == 'SAV'):
    ax.set_xlabel(label, fontsize=fontsize)
  return

mean_df['vpd_crit'] = et_min_vpd(mean_df)

# pft order for panel plots
pft_order = ['DBF', 'EBF', 'MF', 'ENF', 'CSH', 'WSA', 'CRO', 'GRA', 'SAV']

# for xlim on plots with vpd, et max
vpd_max = dfs['95'].vpd.max() + 100.0
vpd_xlim = [0., vpd_max]
t_lims = [dfs['5'].t_a.min() - 1.0, dfs['95'].t_a.max() + 1.0]

def panel_wrapper(_df, function, name, args=()):
  """wrapper that groups by pft and does plot defined by function"""
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*3)
  fig.set_figwidth(fig.get_figwidth()*3)
  for i, pft in enumerate(pft_order):
    print(pft)
    ax = fig.add_subplot(3, 3, i+1)
    function(_df.loc[(_df.pft==pft), :], ax, *args)
  plt.tight_layout()
  plt.savefig('../../doc/paper/%s' % name)
  return
