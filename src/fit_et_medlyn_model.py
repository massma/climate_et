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
from scipy.optimize import least_squares
import codebase.data_io as d_io
import scipy.io as io
import importlib

importlib.reload(pm)
importlib.reload(d_io)

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
  # canopy['g1'] = g_coef[1]
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  data.loc[data.et > 1000., 'et'] = 1000.
  data.loc[data.et <= 0., 'et'] = 0.
  return data['et'] - data['et_obs']

def calc_coef():
  """
  calcualtes best fit coefficients and r2 for each site in ameriflux
  """
  filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])
  _ds = np.ones(len(filenames))*np.nan
  _coef = pd.DataFrame(data={'PFT' : _ds, 'lai' : _ds, 'g1' : _ds,\
                             'count' : _ds, 'r2': _ds}, index=filenames)

  time_start = time.time()
  for filename in filenames[:]:
    print('working on %s' % filename)
    atmos, canopy, data = d_io.load_mat_data(filename)
    _et = data['et_obs']
    if _et.size == 0:
      print('filename %s has no data' % filename)
      continue
    pft = canopy.iloc[0,:].loc['pft']
    try:
      _ = pm.WUE.loc[pft, :]
    except KeyError:
      print("file %s 's pft (%s) has no uWUE, moving on" % (filename, pft))
      continue
    try:
      canopy['g1'] = pm.WUE_MEDLYN.loc[pft, 'g1M']
    except KeyError:
      canopy['g1'] = pm.WUE_MEDLYN.loc[:, 'g1M'].mean()
      print('error, no medlyn coeffieent for %s, pft: %s' % (filename, pft))
      continue

    # result = least_squares(medlyn_fit_et, [1.0, canopy['g1'].iloc[0]],\
    #                        bounds=(0., np.inf),\
    #                        loss='cauchy', args=(atmos, canopy, data))

    result = least_squares(medlyn_fit_et, [1.0],\
                           bounds=(0., 200.),\
                           loss='cauchy', args=(atmos, canopy, data))
    if result['success']:
      _coef.loc[filename, 'lai'] = result['x'][0]
    else:
      print('ERRORR!!!! filename %s failed!' % filename)
      #_coef.loc[filename, 'g1'] = result['x'][1]
      _coef.loc[filename, 'g1'] = canopy['g1'].iloc[0]
    _coef.loc[filename, 'PFT'] = canopy['pft'].iloc[0]
    _coef.loc[filename, 'r2'] = 1. - \
                                np.sum(medlyn_fit_et(result['x'], atmos,\
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
  for key in ['lai', 'g1', 'r2']:
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
coef.to_csv('../dat/site_coef_mm_s_W_m2_medlyn_lai.csv')
statistics = coef.groupby('PFT').apply(generate_coef_stats)
statistics.index = statistics.index.droplevel(1)
outdir = '../dat/adam_mm_s_W_m2_medlyn_lai.csv'
statistics.to_csv(outdir)
#   return coef

# if str(__name__) == '__main__':
#   main()
