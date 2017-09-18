#! ~/edward/bin/python
"""
This script makes some simple analytical plots of vars that
vary in the dET/D_s term
"""
import glob
import os
import pandas as pd
import numpy as np
import metcalcs as met
import seaborn as sns
import resource
from scipy.stats import spearmanr
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import matplotlib as mpl
import codebase.penman_monteith as pm
import util
from datetime import datetime
mpl.rcParams.update(mpl.rcParamsDefault)

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])

df['g_a'] = 1./df['r_a']

time = pd.DatetimeIndex(df.time)
df['hour'] = time.hour
df['jd'] = time.dayofyear
df['year'] = time.year

def get_uwue(_df):
  """lazy way to get uwue"""
  return _df.loc[:, ['uwue']].iloc[0]

def d_et(_df):
  """calcs the full d ET/d Ds, confirmed correct vs df"""
  return _df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

FIG = plt.figure()

def fit_curves(_df, ax):
  """fits curves"""
  jd = _df.jd.values
  lai = _df.lai.values
  def seasonal_lai(x):
    """
    model to fit to lai,
    x: x[0] * sin(jd/365*2.pi + x[1]) + x[2]
    """
    return x[0]*np.sin(jd/365.*2.*np.pi + x[1]) + x[2] - lai
  x0 = [lai.std(), 0., lai.mean()]
  bounds = ([0., 0., 0.], [10., 2.*np.pi, 10.])
  result = least_squares(seasonal_lai, x0, bounds=bounds)
  if result['success']:
    seasonal = seasonal_lai(result['x']) + lai
    df_out = pd.DataFrame(data={'modelled_cycle': seasonal}, index=_df.index)
    for i, _x in enumerate(result['x']):
      df_out['x%d' % i] = _x
    ax.plot(_df.time.values, seasonal, 'k-')
    return df_out
  else:
    print('error, no convergent solution for %s at year %d'\
          % (_df.site.iloc[0], int(_df.year.iloc[0])))
    return np.ones(3)*np.nan
  return

def plot_lai(_df, fit_plot=False, suff=''):
  ax = FIG.add_subplot(111)
  ax.scatter(_df.time.values, _df.lai.values, s=1)
  FIG.autofmt_xdate()
  x = _df.groupby('year').apply(fit_curves, ax)
  util.test_savefig('%s/climate_et/lai_plots/%s_%s%s.png'\
                    % (os.environ['PLOTS'], _df.pft.iloc[0],\
                       _df.site.iloc[0], suff))
  FIG.clf()
  return x

# # below is to generate datta
# x = df.groupby('site').apply(plot_lai, fit_plot=True)
# full_df = pd.concat([df, x], axis=1)
# full_df.to_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
#                   % os.environ['DATA'])

df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA'])

bestfit = df.loc[df.x0.argmax()]

# well-fit subset to test with edward approach
subset = df.loc[(df.year == 2003) & (df.site == 'US-MMS'), :]
plot_lai(subset, suff='_subset')

