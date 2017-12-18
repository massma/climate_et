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
import metcalcs as met
import codebase.data_calc as d_calc
import pickle

mpl.use('Pdf')
mpl.rcParams.update(mpl.rcParamsDefault)
importlib.reload(util)
importlib.reload(plot_tools)
importlib.reload(d_calc)

# grab sites actually used in analysis
with open('%s/changjie/diagnosed_data.pkl' % os.environ['DATA'],\
          mode='rb') as file:
  dfs = pickle.load(file)
df = dfs['full']
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

# df.index = df.reset_index(drop=True)
