#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import glob
import time
import os
import importlib
import pandas as pd
import codebase.penman_monteith as pm
import codebase.data_io as d_io
import numpy as np

importlib.reload(d_io)
importlib.reload(pm)

SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

H2O = 18.01528e-3 #molecular mass kg/mol
def leaf_vpd(atmos, canopy, lai):
  """calculates the leaf term in dET/dVPD (see doc folder)"""
  return -atmos['gamma']*atmos['c_a']*pm.LV/\
    (lai*1.6*pm.R_STAR*canopy['uwue'])\
    *(2.*canopy['g1'] + np.sqrt(atmos['vpd']))\
    /(2.*(canopy['g1'] + np.sqrt(atmos['vpd']))**2)

def calc_lai(atmos, canopy, et_obs):
  """cals lai given structure, with an explcit et_obs"""
  return - atmos['gamma']*atmos['c_a']*pm.LV\
                *np.sqrt(atmos['vpd'])*atmos['p_a']\
                /(atmos['r_a']\
                  *(et_obs*(atmos['gamma']+atmos['delta'])\
                    -atmos['delta']*(atmos['r_n']-canopy['g_flux'])\
                    -1./atmos['r_a']*atmos['rho_a']*pm.CP*atmos['vpd'])
                  *1.6*pm.R_STAR*(273.15 + atmos['t_a'])\
                  *canopy['uwue']*(1.+canopy['g1']/np.sqrt(atmos['vpd'])))

def et_to_gpp(atmos, canopy):
  """converts et to gpp using WUE"""
  return canopy['uwue']/(np.sqrt(atmos['vpd'])*pm.LV)\
                *pm.penman_monteith_uwue(atmos, canopy)

def d_gpp_d_vpd(atmos, canopy):
  """calcualtes analytical derivative d gpp"""
  out = 1./(2.*pm.LV*atmos['vpd']**1.5*(atmos['delta']+atmos['gamma']))\
        *(-canopy['uwue']*atmos['delta']*(atmos['r_n']-canopy['g_flux'])\
          +1./atmos['r_a']*canopy['uwue']*atmos['rho_a']*pm.CP*atmos['vpd']\
          -canopy['g1']*atmos['gamma']*atmos['c_a']*pm.LV*atmos['p_a']\
          /(atmos['r_a']*canopy['lai']*pm.R_STAR*(273.15+atmos['t_a'])*1.6\
            *(1.+canopy['g1']/np.sqrt(atmos['vpd']))**2))
  return out


def calc_derivative(atmos, canopy, data):
  """adds various derivative fields to data, given atmos and canopy"""
  #calculate LAIs
  canopy['lai'] = calc_lai(atmos, canopy, data['et_obs'])
  data['et_gpp'] = data['gpp_obs']*np.sqrt(atmos['vpd'])/canopy['uwue']
  canopy['lai_gpp'] = calc_lai(atmos, canopy, data['et_gpp'])

  # Now do ET terms
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  data['scaling'] = atmos['p_a']/(atmos['r_a']*(273.15+atmos['t_a'])\
                                  *(atmos['delta'] + atmos['gamma']))
  data['vpd_atm'] = pm.CP/atmos['r_moist']
  data['vpd_leaf'] = leaf_vpd(atmos, canopy, canopy['lai'])
  atmos['vpd'] = atmos['vpd'] + 1.0
  data['et_all'] = pm.penman_monteith_uwue(atmos, canopy)
  atmos['vpd'] = atmos['vpd'] - 1.0
  data['d_et'] = data['scaling']*(data['vpd_leaf'] + data['vpd_atm'])
  data['d_et_vpd_std'] = atmos.vpd.std()*data.d_et # units: W/m2
  data['d_et_vpd_std_leaf'] = atmos.vpd.std()*data.vpd_leaf*data.scaling
  data['d_et_vpd_std_atm'] = atmos.vpd.std()*data.vpd_atm*data.scaling

  #Calculate WUE terms
  data['wue_obs'] = data['gpp_obs']/(data['et_obs']/(pm.LV*H2O)*1.e6)
  data['wue'] = canopy['uwue']/np.sqrt(atmos['vpd'])*H2O/1.e6
  data['d_wue'] = -0.5*canopy['uwue']/atmos['vpd']**1.5*H2O/1.e6
  data['d_wue_vpd_std'] = atmos.vpd.std()*data['d_wue']

  # Now GPP terms
  lai_true = canopy['lai'].copy()
  canopy['lai'] = canopy['lai_gpp'].copy()
  data['gpp'] = et_to_gpp(atmos, canopy)
  atmos['vpd'] = atmos['vpd'] + 1.0
  data['gpp_all'] = et_to_gpp(atmos, canopy)
  atmos['vpd'] = atmos['vpd'] - 1.0
  data['d_gpp'] = d_gpp_d_vpd(atmos, canopy)
  data['d_gpp_vpd_std'] = atmos.vpd.std()*data['d_gpp']

  #retun lai to that used for ET
  canopy['lai'] = lai_true
  return atmos, canopy, data

#def main():
"""wrapper for main script"""


start = time.time()
outdir = '%s/changjie/pandas_data_lai/' % os.environ['DATA']

filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])

time_start = time.time()
for filename in filenames[:]:
  print('working on %s' % filename)
  atmos, canopy, data = d_io.load_mat_data(filename)
  if (data.et_obs.count() > 0) & (canopy.dropna().uwue.count() > 0):
    atmos, canopy, data = calc_derivative(atmos, canopy, data)
    dfout = pd.concat([atmos, canopy, data], axis=1)
    fname = ''.join(filename.split('/')[-1].split('.')[:-1])
    dfout.to_pickle('%s/%s.pkl' % (outdir, fname))
  else:
    print('filename %s is invalid, et count %d and canopy count %d'\
         % (filename, data.et_obs.count(), canopy.dropna().uwue.count()))

print('time was %f s' % ((time.time()-start)))
#return

# if str(__name__) == '__main__':
#   main()
reload_data = True
if reload_data:
  concat_dfs(folder='pandas_data_lai', fname='full_pandas_lai')
  df = pd.read_pickle('%s/changjie/full_pandas_lai.pkl' % os.environ['DATA'])
  meta = {}
  meta['folder_label'] = 'site'
  meta['folder'] = 'hist_plots'
  meta['var'] = 'lai_gpp'
  print(df.shape)
  df = df.groupby('site').apply(site_clean)
  print(df.shape)
  df = clean_df(df)
  df = clean_df(df, var='lai_gpp')
  # test = df.groupby('site').apply(site_clean, 'lai_gpp')
  # test = clean_df(test, var='lai_gpp')
  df.to_pickle('%s/changjie/full_pandas_lai_clean.pkl' % os.environ['DATA'])
  print(df.shape)
  #df.groupby('site').apply(histogram, meta)
  # histogram(df, meta)
  # meta['var'] = 'lai'
  # histogram(df, meta)
