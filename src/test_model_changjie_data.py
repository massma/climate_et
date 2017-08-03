#! ~/edward/bin/python
"""
This module re-does the medlyn fit, allowing GPP modify g0
"""
import time
import os
import glob
import importlib
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io as io

importlib.reload(pm)
WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
           comment='#', delimiter=',')
WUE.index = WUE.PFT
LV = 2.5e6

# def medlyn_fit(_g, *data):
#   """function to fit using leastsq"""
#   _vpd, _gpp, _g_s, pft = data
#   return _g[0]*_gpp*(1. + _g[1]/np.sqrt(_vpd)) - _g_s

def medlyn_fit(_g, *data):
  """function to fit using leastsq"""
  _vpd, _gpp, _g_s, _pft, _et = data
  wue = WUE.loc[_pft, 'u_wue_yearly']*1.e6/12.011
  return _g[0]*wue*_et/LV/np.sqrt(_vpd*10.)\
    *(1. + _g[1]/np.sqrt(_vpd)) - _g_s

#def calc_coef():
"""
calcualtes best fit coefficients and r2 for each site in ameriflux
"""
sitelist = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
sitelist.index = sitelist.Site
filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])
_ds = np.ones(len(filenames))*np.nan
_coef = pd.DataFrame(data={'PFT' : _ds, 'g0' : _ds, 'g1' : _ds,\
                           'count' : _ds, 'r2': _ds}, index=filenames)

time_start = time.time()
for i,filename in enumerate(filenames[:]):
  print('working on %d' % i)
  data = io.loadmat(filename)
  pft = str(np.squeeze(data['cover_type']))
  try:
    _ = WUE.loc[pft, :]
  except KeyError:
    print("file %s 's pft has no uWUE, moving on" % filename)
    continue

  
  atmos = {}
  atmos['t_a'] = np.squeeze(data['TA']) #C
  atmos['rh'] = np.squeeze(data['RH']*100.) #percent
  atmos['h'] = np.squeeze(data['H'])
  atmos['r_n'] = np.squeeze(data['NETRAD'])
  atmos['ustar'] = np.squeeze(data['USTAR']) #m/s
  atmos['p_a'] = np.squeeze(data['PA']*1000.) #Pa
  atmos['u_z'] = np.squeeze(data['WS'])
  
  canopy = {}
  canopy['g_flux'] = np.squeeze(data['G'])
  canopy['height'] = float(sitelist.loc[data['sitecode'],'Canopy_h'])
  canopy['zmeas'] = float(sitelist.loc[data['sitecode'],'Measure_h'])
  canopy['r_s'] = np.squeeze(data['Rs']) #s/m


  et = pm.penman_monteith(atmos, canopy)
  if np.isnan(np.nanmax(et)):
    print('error no et for %s' % filename)
  else:
    plt.figure()
    sns.regplot(x=data['LE'], y=et)
    plt.savefig('%s/%05d_garb.png' % (os.environ['PLOTS'], i))
    plt.show(block=False)

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


