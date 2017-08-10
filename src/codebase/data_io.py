#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import time
import os
import glob
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
from scipy.optimize import leastsq
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
  """
  data = io.loadmat(filename)

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
  canopy['height'] = float(SITELIST.loc[data['sitecode'], 'Canopy_h'])
  canopy['zmeas'] = float(SITELIST.loc[data['sitecode'], 'Measure_h'])
  canopy['r_s'] = np.squeeze(data['Rs']) #s/m

  atmos, canopy = pm.penman_monteith_prep(atmos, canopy)
  canopy.pop('r_s')
  atmos = pd.DataFrame(data=atmos, index=np.squeeze(data['time']))
  canopy = pd.DataFrame(data=canopy, index=np.squeeze(data['time']))

  canopyclmns = canopy.columns
  atmosclmns = atmos.columns
  df = pd.concat([atmos, canopy], axis=1).dropna()
  atmos = df[atmosclmns]
  canopy = df[canopyclmns]
  canopy['pft'] = str(np.squeeze(data['cover_type']))
  return atmos, canopy, data

