#! ~/edward/bin/python
"""
This module does some tests to make sure I didn't do any errors
"""
import os
import time
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from codebase.data_calc import *

plt.close('all')

PLOTDIR = '%s/climate_et/test_plots' % os.environ['PLOTS']

def max_diff(quant1, quant2):
  """caclualte the maximuma boslute difference"""
  return np.nanmean(np.absolute(quant1-quant2))

def median_diff(quant1, quant2):
  """caclualte the maximuma boslute difference"""
  return np.nanmedian(np.absolute(quant1-quant2))

def test_et_model(_df):
  """makes a histogram of calc'd uwue as compated to zhou uwue"""
  plt.figure()
  sns.distplot(_df['sigma'])
  plt.xlabel("sigma") #"(calc'ed uwue)/(zhou's uwue)")
  plt.title("PDF for middle 90\% of data")
  plt.savefig('%s/sigma_hist.png' % PLOTDIR)
  return

def test_et_models(_df):
  """calcs histograms of c.v for all models"""
  plt.close('all')
  for name in ['uwue', 'uwue_bb', 'iwue_bb', 'iwue']:
    plt.figure()
    sns.distplot(_df[name]/_df[name].mean())
    plt.xlabel("cv - %s" % name) #"(calc'ed uwue)/(zhou's uwue)")
    plt.title("PDF - %s" % name)
    plt.savefig('%s/cv_hist_%s.png' % (PLOTDIR, name))
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

def plot_et_curve(mean_row, min_row, max_row, index, _df):
  """plots the et curve given a row of df"""
  # _df = _df.loc[(_df.pft == index)]
  # start = time.time()
  # plt.figure()
  # sns.jointplot(_df.vpd, _df.et_obs, kind='hex')
  # plt.savefig('%s/joint_plt_%s.png' % (PLOTDIR, index))
  # print('sns plot time was %f s' % (time.time()-start))
  vpd = np.linspace(min_row.vpd, max_row.vpd)
  et = pm_et(mean_row, vpd=vpd)
  et_orig = pm_et_orig(mean_row, vpd=vpd)
  # et_iwue = pm_et(mean_row, vpd=vpd, uwue=mean_row.iwue, n=1.0)
  # et_iwue_bb = pm_et(mean_row, vpd=vpd, uwue=mean_row.iwue_bb,\
  #                 n=1.0, kernel=bb_kernel)
  # et_bb = pm_et(mean_row, vpd=vpd, uwue=mean_row.uwue_bb,\
  #                 n=0.5, kernel=bb_kernel)
  plt.figure()
  # start = time.time()
  # _df = _df.loc[((_df.r_net > mean_row.r_net-2.0) &\
  #                (_df.r_net < mean_row.r_net+2.0))]
  # end = time.time()
  # print('for pft %s, time to subset was %f s, and nobs was %d'\
  #     % (index, (end-start), _df.shape[0]))
  # plt.scatter(_df.vpd, _df.et_obs, s=4.0, c='m')
  plt.plot(vpd, et, label='uWUE PM')
  plt.plot(vpd, et_orig, label='Mean GPP')
  # plt.plot(vpd, et_iwue, label='iWUE PM')
  # plt.plot(vpd, et_bb, label='uWUE BB')
  # plt.plot(vpd, et_iwue_bb, label='iWUE BB')
  plt.legend(loc='best')
  plt.title(index)
  plt.xlabel('vpd (Pa')
  plt.ylabel('et (w/m2)')
  plt.savefig('%s/%s_et_vpd_curve_with_pm.png' % (PLOTDIR, index))
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
  cp_df.loc[:, 'iwue'] = mean_df.iwue.loc[_df.pft.iloc[0]]
  cp_df.loc[:, 'iwue_bb'] = mean_df.iwue_bb.loc[_df.pft.iloc[0]]
  cp_df.loc[:, 'uwue_bb'] = mean_df.uwue_bb.loc[_df.pft.iloc[0]]
  models = {}
  models['uwue'] = pm_et(cp_df)
  models['uwue_bb'] = pm_et(cp_df, kernel=bb_kernel, uwue=cp_df.uwue_bb)
  models['iwue'] = pm_et(cp_df, uwue=cp_df.iwue, n=1.0)
  models['iwue_bb'] = pm_et(cp_df, kernel=bb_kernel, uwue=cp_df.iwue_bb, n=1.0)
  models['original'] = pm_et_orig(cp_df,\
                                  lai=mean_df.lai_pm.loc[_df.pft.iloc[0]]) #1.0
  # models['gppfixed'] = _df['et_gppfixed']
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
  bias = {'uwue' : [], 'original' : [], 'gppfixed' : []}
  rmse = {'uwue' : [], 'original' : [], 'gppfixed' :[]}
  bias = {'original' : [], 'uwue' : [], 'uwue_bb' : [],\
          'iwue' : [], 'iwue_bb' : []}
  rmse = copy.deepcopy(bias)

  for pft in pfts:
    _bias, _rmse = compare_et(dfs['full'].loc[(dfs['full'].pft == pft), :],\
                              dfs['mean'])
    for key in bias:
      bias[key].append(_bias[key])
      rmse[key].append(_rmse[key])
  # dfs['full'].groupbymake('pft').apply(compare_et, dfs['mean'])
  fig = plt.figure()
  ax1 = fig.add_subplot(211)
  ax2 = fig.add_subplot(212)
  x = np.arange(len(pfts))
  for key in bias:
    ax1.scatter(x, np.absolute(bias[key]), label=key)
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
                 dfs['95'].loc[index], index, _df)\
   for index in dfs['mean'].index]
  compare_et_wrapper(dfs)
  test_et_models(_df)
  return
