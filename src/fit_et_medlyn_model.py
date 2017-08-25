#! ~/edward/bin/python
"""
This module re-does the medlyn fit, allowing GPP modify g0
"""
import time
import os
import glob
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
from scipy.optimize import leastsq
import codebase.data_io as d_io
import scipy.io as io
import importlib

importlib.reload(pm)
importlib.reload(d_io)

WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
           comment='#', delimiter=',')
WUE.index = WUE.PFT
LV = 2.5e6

WUE_MEDLYN = pd.read_csv('../dat/franks_et_al_table2.csv',\
                         comment='#', delimiter=',')
WUE_MEDLYN.index = WUE_MEDLYN.PFT


SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

def medlyn_fit_et(g_coef, *args):
  """
  wrapper to solve for g using obs et
  """
  atmos, canopy, data = args
  canopy['lai'] = g_coef[0]
  canopy['g1'] = g_coef[1]
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  # data.loc[data.et > 1000., 'et'] = 1000.
  # data.loc[data.et <= 0., 'et'] = 0.
  return data['et'] - data['et_obs']

def calc_coef():
  """
  calcualtes best fit coefficients and r2 for each site in ameriflux
  """
  filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])
  _ds = np.ones(len(filenames))*np.nan
  _coef = pd.DataFrame(data={'PFT' : _ds, 'g0' : _ds, 'g1' : _ds,\
                             'count' : _ds, 'r2': _ds}, index=filenames)

  time_start = time.time()
  for filename in filenames[:]:

    atmos, canopy, data = d_io.load_mat_data(filename)
    _et = data['et_obs']
    if _et.size == 0:
      print('filename %s has no data' % filename)
      continue
    pft = canopy.iloc[0,:].loc['pft']
    try:
      _ = WUE.loc[pft, :]
    except KeyError:
      print("file %s 's pft has no uWUE, moving on" % filename)
      continue
    try:
      g_1 = pm.WUE_MEDLYN.loc[pft, 'g1M']
    except KeyError:
      g_1 = pm.WUE_MEDLYN.loc[:, 'g1M'].mean()

    _g, ier = leastsq(medlyn_fit_et, [1.0, g_1],\
               args=(atmos, canopy, data))
    if (ier <= 4) & (ier > 0):
      _coef.loc[filename, 'g0'] = _g[0]
      _coef.loc[filename, 'g1'] = _g[1]
    _coef.loc[filename, 'PFT'] = canopy['pft'].iloc[0]
    _coef.loc[filename, 'r2'] = 1. - \
                                np.sum(medlyn_fit_et(_g, atmos,\
                                                         canopy, data)**2)\
                                                /np.sum((_et - _et.mean())**2)
    _coef.loc[filename, 'count'] = _et.size
  print('wall time was %f s' % (time.time()-time_start))
  print('r2: %f' % _coef.r2.mean())
  return _coef.dropna()


def generate_coef_stats(_coef):
  """
  calcualted mean and std deviation of each PFT,\
  note should prob group by pft and run lsq
  """
  _dout = {}
  for key in ['g0', 'g1', 'r2']:
    _mean = (_coef[key]*_coef['count']).sum()\
        /_coef['count'].sum()
    _dout['%s_mean' % key] = [_mean]
    _dout['%s_std' % key] = [np.sqrt(((_coef[key] - _mean)**2\
                                      *_coef['count']).sum()\
                                     /_coef['count'].sum())]
  return pd.DataFrame(data=_dout)

# def main():
"""wrapper for main script"""
coef = calc_coef()
coef.to_csv('../dat/site_coef_mm_s_W_m2_medlyn.csv')
statistics = coef.groupby('PFT').apply(generate_coef_stats)
statistics.index = statistics.index.droplevel(1)
outdir = '../dat/adam_mm_s_W_m2_medlyn.csv'
statistics.to_csv(outdir)
#   return coef

# if str(__name__) == '__main__':
#   main()
