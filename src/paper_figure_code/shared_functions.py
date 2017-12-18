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

def et_min_vpd(_df, lai):
  """calculates theoretical max vpd as functoin of -df and lai"""
  c3 = pm.CP/_df.r_moist
  # c1 = n_df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
  c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
  c2 = _df.g1
  sqrt_vpd = (c1 + np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3)
  try:
    sqrt_vpd[sqrt_vpd < 0.] = np.nan
  except TypeError:
    if sqrt_vpd < 0.:
      sqrt_vpd = np.nan
      print("setting to nan")
  return sqrt_vpd**2

def get_pft(_df):
  return _df['pft'].iloc[0]

def term_2(_df, lai, vpd):
  """calculates term 2, takes in a _df of one pft"""
  atmos = {'gamma' : _df.gamma.mean(), 'c_a' : _df.c_a.mean(),\
           'vpd' : vpd}
  if _df.uwue.std() > 1.e-8:
    print('error, uWUE is variable: %f!!!!' % _df.uwue.std())
  elif _df.g1.std() > 1.e-8:
    print('error, g1 is variabile: %f!!!!!' % _df.g1.std())
  canopy = {'uwue' : _df.uwue.mean(), 'g1' : _df.g1.mean()}
  return pm.CP/_df.r_moist.mean() +calc.leaf_vpd(atmos, canopy, lai)

# df.index = df.reset_index(drop=True)
