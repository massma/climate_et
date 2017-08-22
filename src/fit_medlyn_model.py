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


WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
           comment='#', delimiter=',')
WUE.index = WUE.PFT
LV = 2.5e6

SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

# def medlyn_fit(_g, *data):
#   """function to fit using leastsq"""
#   _vpd, _gpp, _g_s, pft = data
#   return _g[0]*_gpp*(1. + _g[1]/np.sqrt(_vpd)) - _g_s

def medlyn_fit_div_et(_g, *data):
  """
  function to fit using leastsq,
  vpd should be in kPa, _g and _g_s in mm/s, and et in q/m2,
  same as medlyn_fit except conductance is divided by ET
  """
  _vpd, _g_s, _pft, _et = data
  wue = WUE.loc[_pft, 'u_wue_yearly']*1.e6/12.011
  return _g[0]*wue/LV/np.sqrt(_vpd*10.)\
    *(1. + _g[1]/np.sqrt(_vpd)) - _g_s

def medlyn_fit(_g, *data):
  """
  function to fit using leastsq,
  vpd should be in kPa, _g and _g_s in mm/s, and et in q/m2
  """
  _vpd, _g_s, _pft, _et = data
  wue = WUE.loc[_pft, 'u_wue_yearly']*1.e6/12.011
  return _g[0]*wue*_et/LV/np.sqrt(_vpd*10.)\
    *(1. + _g[1]/np.sqrt(_vpd)) - _g_s

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

    canopy['r_s'] = data['r_s']
    pft = canopy.iloc[0,:].loc['pft']
    try:
      _ = WUE.loc[pft, :]
    except KeyError:
      print("file %s 's pft has no uWUE, moving on" % filename)
      continue

    vpd = atmos['vpd']/1000.
    # note below is in mol/m2/s
    #g_s = data['Gs']
    # below is mm/s ?
    # g_s = data['g_s']
    #below is mm/s/W/m2
    g_s = data['g_s']/_et
    _g, ier = leastsq(medlyn_fit_div_et, [0.04/_et.mean(), 0.7/_et.mean()],\
               args=(vpd, g_s, pft, _et))
    print(g_s.mean())
    if (ier <= 4) & (ier > 0):
      _coef.loc[filename, 'g0'] = _g[0]
      _coef.loc[filename, 'g1'] = _g[1]
    _coef.loc[filename, 'PFT'] = canopy['pft'].iloc[0]
    _coef.loc[filename, 'r2'] = 1. - \
                                np.sum(medlyn_fit_div_et(_g, vpd,\
                                                  g_s, pft, _et)**2)\
                                                /np.sum((g_s - g_s.mean())**2)
    _coef.loc[filename, 'count'] = _et.size
  print('wall time was %f s' % (time.time()-time_start))
  print('r2: %f' % _coef.r2.mean())
  return _coef.dropna()


def generate_coef_stats(_coef):
  """
  calcualted mean and std deviation of each PFT
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

def main():
"""wrapper for main script"""
coef = calc_coef()
coef.to_csv('../dat/site_coef_mm_s_W_m2_medlyn.csv')
statistics = coef.groupby('PFT').apply(generate_coef_stats)
statistics.index = statistics.index.droplevel(1)
outdir = '../dat/adam_mm_s_W_m2_medlyn.csv'
statistics.to_csv(outdir)
  return coef

# if str(__name__) == '__main__':
#   main()
