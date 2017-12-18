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

def plot_et_curve(mean_row, min_row, max_row, index, _function=pm_et):
  """plots the et curve given a row of df"""
  vpd = np.linspace(min_row.vpd, max_row.vpd)
  et = pm_et(mean_row, vpd=vpd)
  et_orig = pm_et_orig(mean_row, vpd=vpd)
  plt.figure()
  plt.plot(vpd, et, label='New PM')
  plt.plot(vpd, et_orig, label='Mean GPP')
  plt.legend(loc='best')
  plt.title(index)
  plt.xlabel('vpd (Pa')
  plt.ylabel('et (w/m2)')
  plt.savefig('%s/%s_et_vpd_curve.png' % (PLOTDIR, index))
  return


def get_bias(models, _df):
  """calcs bias for all three"""
  bias = {}
  for key in models:
    bias[key] = np.mean(models[key]-_df['et_obs'])
  return bias

def get_rmse(models, _df):
  """calcs rmse for all three"""
  rmse = {}
  for key in models:
    rmse[key] = np.sqrt(np.mean((models[key]-_df['et_obs'])**2))
  return rmse

def compare_et(_df, mean_df):
  """
  prints some stats comparing PM to our new PM,\
  to be fair should probably fit LAI
  """
  print('\n PFT: %s' % _df.pft.iloc[0])
  cp_df = _df.copy()
  cp_df.loc[:, 'uwue'] = mean_df.uwue.loc[_df.pft.iloc[0]]
  cp_df['lai_pm'] = cp_df.lai_pm.mean() # 1.0
  models = {}
  models['new'] = pm_et(cp_df)
  models['original'] = pm_et_orig(cp_df)
  models['pet'] = cp_df['pet']
  bias = get_bias(models, _df)
  rmse = get_rmse(models, _df)
  for key in bias:
    print('%s model bias: %f W/m2' % (key, bias[key]))
  for key in rmse:
    print('%s model rmse: %f W/m2' % (key, rmse[key]))
  #return
  return bias, rmse

def compare_et_wrapper(dfs):
  """organizes figure and wraps compare_et, which compares new to old PM"""
  pfts = list(dfs['mean'].index)
  bias = {'new' : [], 'original' : []}
  rmse = {'new' : [], 'original' : []}
  for pft in pfts:
    _bias, _rmse = compare_et(dfs['full'].loc[(dfs['full'].pft == pft), :],\
                            dfs['mean'])
    for key in bias:
      bias[key].append(_bias[key])
      rmse[key].append(_rmse[key])
  # dfs['full'].groupby('pft').apply(compare_et, dfs['mean'])
  fig = plt.figure()
  ax1 = fig.add_subplot(211)
  ax2 = fig.add_subplot(212)
  x = np.arange(len(pfts))
  for key in bias:
    ax1.scatter(x, bias[key], label=key)
    ax2.scatter(x, rmse[key], label=key)
  ax1.set_xticks([])
  ax1.set_ylabel('Bias (W/m2)')
  ax2.set_xticks(x)
  ax2.set_xticklabels(pfts)
  ax2.set_ylabel('RMSE (W/m2)')
  plt.legend(loc='best', frameon=True)
  plt.savefig('%s/rmse_bias.png' % PLOTDIR)
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
  compare_et_wrapper(dfs)
  return
