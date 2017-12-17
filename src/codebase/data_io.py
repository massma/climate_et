#! ~/edward/bin/python

"""
This module loads in data prepared by Changjie,
and also loads in some data from tables from Zhou and Franks et. al.
"""
import os
import glob
import pandas as pd
import numpy as np
import scipy.io as io

LV = 2.5e6
WUE_MEDLYN = pd.read_csv('../dat/franks_et_al_table2.csv',\
                         comment='#', delimiter=',')
WUE_MEDLYN.index = WUE_MEDLYN.PFT
# convert from sqrt(kPa) to sqrt(Pa)
WUE_MEDLYN.loc[:, 'g1M'] = WUE_MEDLYN.loc[:, 'g1M']*np.sqrt(1000.)

WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
                  comment='#', delimiter=',')
WUE.index = WUE.PFT
# convert from g C to micromol (units of c_s),
# and from sqrt(hPa) to sqrt(PA)
# and from kg H20 to joules
WUE.loc[:, 'u_wue_yearly'] = WUE.loc[:, 'u_wue_yearly']\
                             *1.e6/12.011*np.sqrt(100.)/LV


SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

def add_atmos_dict(data_out, data):
  """generates an atmos dictionary"""
  data_out['t_a'] = np.squeeze(data['TA']) #C
  data_out['rh'] = np.squeeze(data['RH']*100.) #percent
  data_out['sensible'] = np.squeeze(data['H']) # w/m2
  data_out['r_n'] = np.squeeze(data['NETRAD']) # w/m2
  data_out['ustar'] = np.squeeze(data['USTAR']) #m/s
  data_out['p_a'] = np.squeeze(data['PA']*1000.) #Pa
  data_out['u_z'] = np.squeeze(data['WS']) # m/s
  data_out['u_z'][data_out['u_z'] <= 0.] = np.nan
  if data['flag_Ca'] == 1:
    data_out['c_a'] = np.squeeze(data['Ca'])
  else:
    data_out['c_a'] = np.ones(data_out['p_a'].shape)*400. #ppm
  return data_out

def add_canopy_dict(data_out, data):
  """generates a canopy dict from loaded data"""
  data_out['g_flux'] = np.squeeze(data['G'])
  data_out['height'] = float(SITELIST.loc[data['sitecode'], 'Canopy_h'])
  data_out['zmeas'] = float(SITELIST.loc[data['sitecode'], 'Measure_h'])
  data_out['pft'] = str(np.squeeze(data['cover_type']))
  return data_out

def gen_time(data):
  """creates time, but not sure if start or end see Read_data.m"""
  datestr = [np.datetime64('%d-%02d-%02dT%02d:00'\
             % (int(year), int(month), int(day), int(hours)))\
             for year, month, day, hours in \
             zip(data['year'], data['month'], data['day'], data['hours'])]
  #time = np.datetime64()
  return datestr

def load_file(filename):
  """loads up a matlab file from changjie"""
  data = io.loadmat(filename)
  data_out = {}
  # below mm/s
  data_out['et_obs'] = np.squeeze(data['LE'])
  data_out['swc'] = np.squeeze(data['SWC'])
  data_out['gpp_obs'] = np.squeeze(data['GEP'])
  data_out['r_a_changie'] = np.squeeze(data['Ra']) # s/m

  data_out = add_atmos_dict(data_out, data)
  data_out = add_canopy_dict(data_out, data)
  data_out['time'] = gen_time(data)
  data_out = pd.DataFrame(data=data_out, index=np.squeeze(data['time']))
  data_out = data_out.dropna()
  data_out.loc[data_out['et_obs'] <= 0.0, 'et_obs'] = np.nan
  try:
    pft = data_out.iloc[0, :].loc['pft']
  except IndexError:
    print('file %s has no data, exiting' % filename)
    pft = 'nan'
    return None
  try:
    data_out['g1'] = WUE_MEDLYN.loc[pft, 'g1M']
  except KeyError:
    data_out['g1'] = np.nan
    print('error, no medlyn coeffieent for %s, pft: %s, setting to none'\
          % (filename, pft))
    return None
  try:
    data_out['uwue_zhou'] = WUE.loc[pft, 'u_wue_yearly']
  except KeyError:
    data_out['uwue_zhou'] = np.nan
  data_out['site'] = ''.join(filename.split('/')[-1].split('.')[:-1])
  data_out = data_out.dropna()
  return data_out

def load_mat_data():
  """loads all mat data provided by changjie, outputing dataframe"""
  filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])
  data_list = [load_file(filename) for filename in filenames]
  data_list = [df for df in data_list if df is not None]
  all_site_data = pd.concat(data_list)
  all_site_data = all_site_data.reset_index()
  return all_site_data
