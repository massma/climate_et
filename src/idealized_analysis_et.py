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
import util
import importlib
mpl.rcParams.update(mpl.rcParamsDefault)

importlib.reload(util)

df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA'])
df['g_a'] = 1./df['r_a']
df['r_net'] = df['r_n'] - df['g_flux']

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

def d_et_d_r_net(_df):
  """derivantve w.r.t net radiation"""
  return _df['delta']/(_df['delta']+_df['gamma'])

def d_et_d_lai(_df):
  """calc derivative d et/ dlai"""
  return _df['g_a']*_df['p_a']*_df['gamma']*_df['c_a']\
    *np.sqrt(_df['vpd'])*pm.LV\
    /(_df['t_a_k']*(_df['delta']+_df['gamma'])*_df['lai']**2\
      *pm.R_STAR*1.6*_df['uwue']*(1. + _df['g1']/np.sqrt(_df['vpd'])))

def d_et_d_g_a(_df):
  """calc derivative w.r.t. g_a"""
  return _df['p_a']/(_df['t_a_k']*(_df['delta'] + _df['gamma']))\
    *(pm.CP*_df['vpd']/_df['r_moist']\
      -_df['gamma']*_df['c_a']*np.sqrt(_df['vpd'])*pm.LV\
      /(_df['lai']*pm.R_STAR*1.6*_df['uwue']\
        *(1. + _df['g1']/np.sqrt(_df['vpd']))))

def d_et_d_delta(_df):
  """calc derivative w.r.t. delta"""
  return (_df['gamma']*(_df['r_n']-_df['g_flux'])\
          -_df['g_a']*_df['p_a']/_df['t_a_k']\
          *(pm.CP*_df['vpd']/_df['r_moist']\
            -_df['gamma']*_df['c_a']*np.sqrt(_df['vpd'])*pm.LV\
            /(_df['lai']*pm.R_STAR*1.6*_df['uwue']\
              *(1. + _df['g1']/np.sqrt(_df['vpd'])))))\
              /(_df['delta'] + _df['gamma'])**2

jacobians = {'vpd' : d_et,\
             'lai' : d_et_d_lai,\
             'seasonal_lai' : d_et_d_lai,\
             'residual_lai' : d_et_d_lai,\
             'g_a' : d_et_d_g_a,\
             'delta' : d_et_d_delta,\
             'r_net' : d_et_d_r_net}

def site_analysis(_df, swc=False):
  """caluclates jacobians for each site"""
  out = {}
  for key in jacobians:
    out[key] = jacobians[key](_df)
    if swc:
      out['swc_%s' % key] = np.polyfit(_df['swc'], _df[key], deg=1)[0]*out[key]

  out = pd.DataFrame(data=out, index=[_df.index])
  return out

def get_pft(_df):
  return _df['pft'].iloc[0]

pft = df.groupby('site').apply(get_pft)
mean = df.groupby('site').mean()
std = df.groupby('site').std()
jacobian = site_analysis(mean)

columns = ['delta', 'g_a', 'seasonal_lai', 'residual_lai',\
           'lai', 'vpd', 'd_et', 'r_net']

def remove_mean(_df):
  """strips mean of column"""
  # print(_df.shape)
  _mean = _df.mean()
  # print(_df.mean().index)
  for column in _mean.index:
    _df[column] = _df[column]-_mean[column]
  return _df

nonlinear = site_analysis(df, swc=False)
deviation = df.loc[:, nonlinear.columns]
deviation['site'] = df['site']
deviation = deviation.groupby('site').apply(remove_mean)
deviation.drop('site', axis=1, inplace=True)

non_vector = np.absolute(nonlinear*deviation)
non_vector['site'] = df['site']

def debug(_df):
  """figure out what the hell is going on"""
  print(_df.shape)
  print(np.mean(_df.loc[:, float_columns]))
  return

non_mean = non_vector.groupby('site').mean()
non_mean['pft'] = pft
non_vector.drop('site', axis=1, inplace=True)
non_vector['pft'] = df['pft']
non_mean_pft = non_vector.groupby('pft').mean()


# std = std.loc[:, columns]

def dot_variability(_jacobian, _std):
  """for multiplying jacovian and std, handles swc columns"""
  out = pd.DataFrame(index=_jacobian.index)
  for name in _jacobian.columns:
    out[name] = _jacobian[name]*_std[name]
  return out

# error = np.absolute(np.sqrt((dot_variability(jacobian**2, std**2))\
#                             .sum(axis=1))-std['d_et'])
# rel_error = error/std['d_et']
# # print(error)
# # print(rel_error)

var_vector = np.absolute(dot_variability(jacobian, std))
var_vector['pft'] = pft

def pft_plot(_df, folder_label=''):
  """acts on df grouped by pft, plots partial derivatives of all vars"""
  pft = _df['pft'].iloc[0]
  print('pft %s, has %d sites' % (pft, _df.vpd.count()))
  # _df.drop('pft')
  fig = plt.figure()
  fig.set_figwidth(fig.get_figwidth()*2)
  ax = fig.add_subplot(111)
  _df.boxplot(ax=ax)
  ax.set_title('pft: %s' % pft)
  util.test_savefig('%s/climate_et/et_jacobian/%s/%s.png'\
                    % (os.environ['PLOTS'], folder_label, pft))
  plt.close('all')
  return

var_vector.groupby('pft').apply(pft_plot)
non_mean.groupby('pft').apply(pft_plot, 'non_linear')

def site_plot(_df, folder_label=''):
  """acts on df grouped by pft, plots partial derivatives of all vars"""
  try:
    pft = str(_df['pft'])
    _ds = pd.Series(_df.drop('pft'))
  except KeyError:
    print('key error')
    pft = _df.name
    _ds = pd.Series(_df)
  site = _df.name

  fig = plt.figure()
  fig.set_figwidth(fig.get_figwidth()*2)
  ax = fig.add_subplot(111)

  _ds.plot(kind='bar',ax=ax)
  ax.set_title('pft: %s' % pft)
  util.test_savefig('%s/climate_et/et_jacobian/site/%s/%s_%s.png'\
                    % (os.environ['PLOTS'], folder_label, pft, site))
  plt.close('all')
  return

non_mean_pft.apply(site_plot, folder_label='mean_non-linear', axis=1)
non_mean.apply(site_plot, folder_label='non-linear', axis=1)
var_vector.apply(site_plot, axis=1)
