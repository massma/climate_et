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
mpl.rcParams.update(mpl.rcParamsDefault)

names = ['r_a', 'p_a', 't_a', 'delta', 'gamma', 'r_moist', 'lai', 'c_a',\
         'uwue', 'g1', 'vpd', 'pft', 'site', 'd_et']
df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA']).loc[:, names]
df['g_a'] = 1./df['r_a']

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

def partial_p_a(_df):
  """calculates partial derivative of d-et/d_Ds w.r.t. P"""
  return _df['g_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_g_a(_df):
  """w.r.t. g_a"""
  return _df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_t_a(_df):
  """w.r.t. t_a"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)**2*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_delta(_df):
  """w.r.t. delta"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta'])**2)*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_gamma(_df):
  """w.r.t. gamma"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta'])**2)*\
    (pm.CP/_df['r_moist']+_df['delta']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_r_moist(_df):
  """w.r.t. r_moist"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']**2)

def partial_c_a(_df):
  """w.r.t. c_a"""
  return _df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (-_df['gamma']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_lai(_df):
  """w.r.t. lai"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']**2*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def partial_vpd(_df):
  """w.r.t. vpd"""
  return -_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (3.*_df['g1']+np.sqrt(_df['vpd']))\
     /(4.*np.sqrt(_df['vpd'])*(_df['g1']+np.sqrt(_df['vpd']))**3))

jacobians = {'p_a' : partial_p_a, 'g_a' : partial_g_a, 't_a' : partial_t_a,\
             'delta': partial_delta, 'gamma' : partial_gamma,\
             'r_moist' : partial_r_moist, 'c_a' : partial_c_a,\
             'lai': partial_lai, 'vpd' : partial_vpd}

def site_analysis(_df):
  """caluclates jacobians for each site"""
  out = {}
  for key in jacobians:
    out[key] = jacobians[key](_df)
  out = pd.DataFrame(data=out, index=[_df.index])
  return out

mean = df.groupby('site').mean()
std = df.groupby('site').std()
jacobian = site_analysis(mean)

_mean = mean.loc[:, jacobian.columns]
_std = std.loc[:, jacobian.columns]

error = np.absolute(np.sqrt((jacobian**2*_std**2).sum(axis=1))-std['d_et'])
rel_error = error/std['d_et']
print(rel_error)


def pft_plot(_df):
  """acts on df grouped by pft, plots partial derivatives of all vars"""
  


for 
