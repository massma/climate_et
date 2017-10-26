#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
these are the functions originally used in ../calc_et_vpd_derivative.py
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



def d_et_d_g_a(_df):
  """calc derivative w.r.t. g_a"""
  return _df['p_a']/(_df['t_a_k']*(_df['delta'] + _df['gamma']))\
    *(pm.CP*_df['vpd']/_df['r_moist']\
      -_df['gamma']*_df['c_a']*np.sqrt(_df['vpd'])*pm.LV\
      /(_df['lai']*pm.R_STAR*1.6*_df['uwue']\
        *(1. + _df['g1']/np.sqrt(_df['vpd']))))

def d_et_d_delta(_df):
  """calc derivative w.r.t. delta"""
  return (_df['gamma']*(_df['r_n']-_df['g_flux'])\
          -_df['g_a']*_df['p_a']/_df['t_a_k']\
          *(pm.CP*_df['vpd']/_df['r_moist']\
            -_df['gamma']*_df['c_a']*np.sqrt(_df['vpd'])*pm.LV\
            /(_df['lai']*pm.R_STAR*1.6*_df['uwue']\
              *(1. + _df['g1']/np.sqrt(_df['vpd'])))))\
              /(_df['delta'] + _df['gamma'])**2


def scaling(atmos):
  """calcualtes the scaling term"""
  return 2.*atmos['p_a']/(atmos['r_a']*(273.15+atmos['t_a'])\
                       *(atmos['delta'] + atmos['gamma']))


def d_et_d_lai(_df):
  """calc derivative d et/ dlai"""
  return _df['g_a']*_df['p_a']*_df['gamma']*_df['c_a']\
    *np.sqrt(_df['vpd'])*pm.LV\
    /(_df['t_a_k']*(_df['delta']+_df['gamma'])*_df['lai']**2\
      *pm.R_STAR*1.6*_df['uwue']*(1. + _df['g1']/np.sqrt(_df['vpd'])))

def calc_derivative(atmos, canopy, data):
  """adds various derivative fields to data, given atmos and canopy"""
  #calculate LAIs
  canopy['lai'] = calc_lai(atmos, canopy, data['et_obs'])
  data['et_gpp'] = data['gpp_obs']*np.sqrt(atmos['vpd'])/canopy['uwue']
  canopy['lai_gpp'] = calc_lai(atmos, canopy, data['et_gpp'])

  # Now do ET terms
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  data['scaling'] = scaling(atmos)
  data['vpd_atm'] = pm.CP/atmos['r_moist']
  data['vpd_leaf'] = leaf_vpd(atmos, canopy, canopy['lai'])
  atmos['vpd'] = atmos['vpd'] + 1.0
  data['et_all'] = pm.penman_monteith_uwue(atmos, canopy)
  atmos['vpd'] = atmos['vpd'] - 1.0
  data['d_et'] = data['scaling']*(data['vpd_leaf'] + data['vpd_atm'])
  data['d_et_vpd_std'] = atmos.vpd.std()*data.d_et # units: W/m2
  data['d_et_vpd_std_leaf'] = atmos.vpd.std()*data.vpd_leaf*data.scaling
  data['d_et_vpd_std_atm'] = atmos.vpd.std()*data.vpd_atm*data.scaling
  concatted = pd.concat([atmos, canopy], axis=1)
  data['d_et_d_lai'] = d_et_d_lai(concatted)
  data['d_et_d_g_a'] = d_et_d_g_a(concatted)
  data['d_et_d_delta'] = d_et_d_delta(concatted)

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

#### below fucntions added from jacobian/varaibility plot
def get_uwue(_df):
  """lazy way to get uwue"""
  return _df.loc[:, ['uwue']].iloc[0]

def d_et(_df):
  """calcs the full d ET/d Ds, confirmed correct vs df"""
  return 2.*_df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

def d_et_d_r_net(_df):
  """derivantve w.r.t net radiation"""
  return _df['delta']/(_df['delta']+_df['gamma'])

