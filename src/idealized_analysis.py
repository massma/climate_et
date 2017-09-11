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
import matplotlib.pyplot as plt
import matplotlib as mpl
import codebase.penman_monteith as pm
mpl.rcParams.update(mpl.rcParamsDefault)

names = ['r_a', 'p_a', 't_a', 'delta', 'gamma', 'r_moist', 'lai', 'c_a',\
         'uwue', 'g1', 'vpd', 'pft', 'site', 'd_et']
df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA']).loc[:, names]
df['g_a'] = 1./df['r_a']

def get_uwue(_df):
  """lazy way to get uwue"""
  return _df.loc[:, ['uwue']].iloc[0]


df['vpd_term'] = (2.*df['g1'] + np.sqrt(df['vpd']))/(\
                 2.*(df['g1'] + np.sqrt(df['vpd']))**2)
df['vpd_scale'] = df['gamma']*df['c_a']*pm.LV/\
                 (1.6*pm.R_STAR*df['uwue'])
df['lai_term'] = 1./df['lai']
df['air_term'] = pm.CP/df['r_moist']
df['scaling'] = df['g_a']*df['p_a']/\
                ((df['t_a']+273.15)*(df['delta'] + df['gamma']))
df['d_et_2'] = df['scaling']*(df['air_term']\
                              -df['vpd_scale']*df['lai_term']*df['vpd_term'])

mean = df.groupby('pft').mean()
std = df.groupby('pft').std()

def variability(_df):
  """calcs variability in a comparable way"""
  out = {}
  _paren_mean = (_df['air_term']\
                 -_df['vpd_scale']*df['lai_term']*_df['vpd_term']).mean()
  _scale_mean = _df['scaling'].mean()
  out['scale_var'] = _df['scaling'].std()*_paren_mean
  out['air_var'] = _scale_mean*_df['air_term'].std()
  out['vpd_term_var'] = _scale_mean*\
                        (_df['vpd_scale']*df['lai_term']).mean()\
                        *_df['vpd_term'].std()
  out['vpd_scale_var'] = _scale_mean*\
                        _df['vpd_scale'].std()*df['lai_term'].mean()\
                        *_df['vpd_term'].mean()
  out['lai_term_var'] = _scale_mean*\
                        _df['vpd_scale'].mean()*df['lai_term'].std()\
                        *_df['vpd_term'].mean()
  out['vpd_var'] = _scale_mean*(_df['vpd_scale']\
                                *df['lai_term']*_df['vpd_term']).std()
  out['pft'] = _df['pft'].iloc[0]
  out = pd.DataFrame(data=out, index=[_df.pft.iloc[0]])
  return out

vari = df.groupby('site').apply(variability)
vari.index = vari.index.droplevel(1)
test = vari.groupby('pft').mean()

print(test)
