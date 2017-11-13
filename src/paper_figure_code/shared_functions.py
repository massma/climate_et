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
import codebase.penman_monteith as pm
import codebase.calc_tools as calc
from sympy import Symbol, sqrt, series, latex, limit


#try some analytis with sympy
g1 = Symbol('g_1')
x = Symbol('\frac{\sqrt{D}}{g_1}')
func = g1*(2 + x)/(2*g1**2*(1 + x)**2)
print(latex(series(func, x, x0=0., dir='+', n=4)))
# x = Symbol('\sqrt{D}')
# func  = (2*g1 + x)/(2*(g1 + x)**2)
# print(latex(series(func, x, x0=0., dir='+')))

mpl.use('Pdf')
mpl.rcParams.update(mpl.rcParamsDefault)
importlib.reload(util)
importlib.reload(plot_tools)
# grab sites actually used in analysis

df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA']).loc[:, 'site'].drop_duplicates()


# get metadata
site_list = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                   'Site_list_(canopy_height).csv' % (os.environ['DATA']))
# subset by sites actually used
site_list.index = site_list.Site
site_list = site_list.loc[df.values]

# below is the label

# below commented df should (and is - I checked) equicalent to fix_scaling_clean
# df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
#                     % os.environ['DATA'])
# # adjustments due to errors found on 10/25
# df['scaling'] *= 2.
# df['d_et'] *= 2.
# df['d_et_vpd_std_leaf'] *= 2.
# df['d_et_vpd_std'] *= 2.
# df['d_et_vpd_std_atm'] *= 2.
# df2 = pd.read_pickle('%s/changjie/full_pandas_fix_scaling_clean.pkl'\
#                      % os.environ['DATA'])
# print('old df shape:', df.shape)
# print('new df shape:', df2.shape)
# _df = df.apply(pd.to_numeric, errors='coerce')
# _df2 = df2.apply(pd.to_numeric, errors='coerce')
# difference = np.nanmax(np.absolute(_df-_df2), axis=0)
# print(difference.shape)
# for column1, column2,  _diff in zip(_df.columns, _df2.columns, difference):
#   print(column1, column2, _diff)

df = pd.read_pickle('%s/changjie/full_pandas_fix_scaling_clean.pkl'\
                    % os.environ['DATA'])

df['g_a'] = 1./df['r_a']
df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
df['d_et_atm'] = df['scaling']*df['vpd_atm']
df['uwue_norm'] = df.uwue/pm.LV
df['r_net'] = df['r_n'] - df['g_flux']
jacobians = {'vpd' : calc.d_et,\
             'lai' : calc.d_et_d_lai,\
             'seasonal_lai' : calc.d_et_d_lai,\
             'residual_lai' : calc.d_et_d_lai,\
             'g_a' : calc.d_et_d_g_a,\
             'delta' : calc.d_et_d_delta,\
             'r_net' : calc.d_et_d_r_net}


def et_min_vpd(_df, lai):
  """calculates theoretical max vpd as functoin of -df and lai"""
  c3 = pm.CP/_df.r_moist
  c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
  c2 = _df.g1
  sqrt_vpd = (c1 + np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3)
  try:
    sqrt_vpd[sqrt_vpd < 0.] = np.nan
  except TypeError:
    if sqrt_vpd < 0.:
      sqrt_vpd = np.nan
  return sqrt_vpd**2
