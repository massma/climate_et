#! ~/edward/bin/python
"""
This module does some tests to make sure I didn't do any errors
"""
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from codebase.data_calc import *

PLOTDIR = '%s/climate_et/test_plots' % os.environ['PLOTS']

def max_diff(quant1, quant2):
  """caclualte the maximuma boslute difference"""
  return np.nanmean(np.absolute(quant1-quant2))

def median_diff(quant1, quant2):
  """caclualte the maximuma boslute difference"""
  return np.nanmedian(np.absolute(quant1-quant2))

def test_et_model(_df):
  """makes a histogram of calc'd uwue as compated to zhou uwue"""
  sigma = _df['uwue']/_df['uwue_zhou']
  plt.figure()
  sns.distplot(sigma)
  plt.xlabel("(calc'ed uwue)/(zhou's uwue)")
  plt.title("Middle 90\% of data")
  plt.savefig('%s/sigma_hist.png' % PLOTDIR)
  return

def test_pm(_df):
  """makes sure that modelled pm matches _df"""
  print('maximum difference for ET model: %f'\
        % (max_diff(_df['et'], _df['et_obs'])))

def numeric_d_et(_df, vpd=None):
  """calcualtes numerical d_et"""
  if vpd is None:
    vpd = _df['vpd']
  return pm_et(_df, vpd=(vpd+1.0))-pm_et(_df, vpd=vpd)

def test_d_et_model(_df):
  """this runs some tests to make sure our derivatives are correct"""
  d_et_num = numeric_d_et(_df)
  print('Max. difference between derivative and numerical',\
        max_diff(_df['d_et']/2.0, d_et_num))
  return

def plot_et_curve(mean_row, min_row, max_row, index):
  """plots the et curve given a row of df"""
  vpd = np.linspace(min_row.vpd, max_row.vpd)
  et = pm_et(mean_row, vpd=vpd)
  plt.figure()
  plt.plot(vpd, et)
  plt.title(index)
  plt.xlabel('vpd (Pa')
  plt.ylabel('et (w/m2)')
  plt.savefig('%s/%s_et_vpd_curve.png' % (PLOTDIR, index))
  return

def run_all_tests(dfs):
  """runs all tests"""
  _df = dfs['full']
  test_d_et_model(_df)
  test_et_model(_df)
  test_pm(_df)
  [plot_et_curve(dfs['mean'].loc[index],\
                 dfs['5'].loc[index],\
                 dfs['95'].loc[index], index) for index in dfs['mean'].index]
  return
