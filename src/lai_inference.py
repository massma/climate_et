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

site = 'FR-Gri'
df = pd.read_pickle('%s/changjie/pandas_data_lai/%s.pkl'\
                    % (os.environ['DATA'], site)) #.loc[:, names]

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
