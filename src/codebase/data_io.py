#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import os
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
import scipy.io as io


SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site


def load_mat_data(filename):
  """
  takes a filename (matlab file from changjie),
  and returns data structures for penman_monteith,
  as well as a dictionary for the mat file
  _data is misc. data to be output
  """
  data = io.loadmat(filename)
  _data = {}
  # below mm/s
  _data['et_obs'] = np.squeeze(data['LE'])
  _data['swc'] = np.squeeze(data['SWC'])

  atmos = {}
  atmos['t_a'] = np.squeeze(data['TA']) #C
  atmos['rh'] = np.squeeze(data['RH']*100.) #percent
  atmos['h'] = np.squeeze(data['H'])
  atmos['r_n'] = np.squeeze(data['NETRAD'])
  atmos['ustar'] = np.squeeze(data['USTAR']) #m/s
  atmos['p_a'] = np.squeeze(data['PA']*1000.) #Pa
  atmos['u_z'] = np.squeeze(data['WS'])
  if data['flag_Ca'] == 1:
    atmos['c_a'] = np.squeeze(data['Ca'])
  else:
    atmos['c_a'] = np.ones(atmos['p_a'].shape)*400. #ppm

  canopy = {}
  canopy['g_flux'] = np.squeeze(data['G'])
  canopy['height'] = float(SITELIST.loc[data['sitecode'], 'Canopy_h'])
  canopy['zmeas'] = float(SITELIST.loc[data['sitecode'], 'Measure_h'])
  canopy['r_s'] = np.squeeze(data['Rs']) #s/m
  atmos, canopy = pm.penman_monteith_prep(atmos, canopy)
  atmos['vpd'][atmos['vpd'] <= 0.] = np.nan
  canopy.pop('r_s')
  atmos = pd.DataFrame(data=atmos, index=np.squeeze(data['time']))
  canopy = pd.DataFrame(data=canopy, index=np.squeeze(data['time']))
  _data = pd.DataFrame(data=_data, index=np.squeeze(data['time']))
  _dataclmns = _data.columns
  canopyclmns = canopy.columns
  atmosclmns = atmos.columns
  df = pd.concat([atmos, canopy, _data], axis=1).dropna()
  _data = df[_dataclmns].copy()
  atmos = df[atmosclmns].copy()
  canopy = df[canopyclmns].copy()
  canopy['pft'] = str(np.squeeze(data['cover_type']))
  # below is s/m
  _data['r_s'] = ((((atmos['delta']*(atmos['r_n']-canopy['g_flux'])+\
                    (1012.*atmos['rho_a']*atmos['vpd'])/atmos['r_a'])\
                   /np.squeeze(_data['et_obs'])-atmos['delta'])\
                  /atmos['gamma'])-1.)*atmos['r_a']
  #below is mm/s
  _data['g_s'] = 1./_data['r_s']*1000.
  _data.loc[_data['g_s'] < 0., 'g_s'] = np.nan
  _data.loc[_data['g_s'] > 100., 'g_s'] = np.nan
  _data.loc[_data['et_obs'] <= 0.0, 'et_obs'] = np.nan
  _data = _data.dropna()
  atmos = atmos.loc[_data.index, :]
  canopy = canopy.loc[_data.index, :]
  return atmos, canopy, _data
