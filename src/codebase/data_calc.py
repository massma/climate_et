#! ~/edward/bin/python

"""
This module calculates all diagnostics used in the paper
"""
import numpy as np

#geophysical constats
VP_FACTOR = 100. #convert hPa -> Pa
K = 0.41 # vonkarmans constant
CP = 1012. # specific heat air, got from Changjie (I usually use 1004)
GAMMA = 66. # psychrometric constant
LV = 2.5e6
R_AIR = .0289645 # mean kg/mol dry air
R_STAR = 8.3144598 #J /mol/K
R_DRY = 287.058
G = 9.81 #gravity constant, m/s2


def uwue(_df):
  """
  calcs uwue from _df with obs
  note possible issue with c_a instead of c_s
  """
  _df['uwue'] = -_df['g_a']*_df['gamma']*_df['c_a']\
                *np.sqrt(_df['vpd'])*_df['p_a']\
                /((_df['et_obs']*(_df['delta']+_df['gamma'])\
                   -_df['delta']*_df['r_net']\
                   -_df['g_a']*_df['rho_a']*CP*_df['vpd'])\
                  *1.6*R_STAR*_df['t_a_k']*(1.0+_df['g1']/np.sqrt(_df['vpd'])))
  return _df

def clean_df(_df):
  """
  cleans df, for now removing lower and upper 5 percent of values,
  based on calc'd uwue. Could look into more sophisticated outlier
  identification (seems like from histogram there are clear ouliers).
  """
  clean_df = _df.loc[((_df.uwue < _df.uwue.quantile(q=0.95)) &\
                      (_df.uwue > _df.uwue.quantile(q=0.05)))]
  return clean_df

def pm_et(_df, vpd=None):
  """calculates et using our new penman-monteith"""
  if vpd is None:
    vpd = _df['vpd']
  _df['et'] = (_df['delta']*_df['r_net']\
               +_df['g_a']*_df['p_a']/_df['t_a_k']\
               *(CP*vpd/_df['r_moist']\
                 -_df['gamma']*_df['c_a']*np.sqrt(vpd)\
                 /(R_STAR*1.6*_df['uwue']\
                   *(1.0+_df['g1']/np.sqrt(vpd)))))\
                   /(_df['delta']+_df['gamma'])
  return _df

def sign(_df, vpd=None):
  """calculates the 'sign' term"""
  if vpd is None:
    vpd = _df['vpd']
  _df['sign'] = CP/_df['r_moist']\
                -_df['gamma']*_df['c_a']\
                /(1.6*R_STAR*_df['uwue'])\
                *((2.*_df['g1']+np.sqrt(vpd))\
                  /(2.0*(_df['g1']+np.sqrt(vpd))**2))
  return _df

def scaling(_df):
  """
  calculates the 'scaling' term,
  note 2.0 factor is from esat-rh decomposition,
  but need to double check this
  """
  _df['scaling'] = 2.0*_df['g_a']*_df['p_a']\
                   /(_df['t_a_k']*(_df['delta'] + _df['gamma']))
  return _df

def d_et_d_vpd(_df, vpd=None):
  """calculates the derivative of et w.r.t. vpd"""
  if vpd is None:
    vpd = _df['vpd']
  _df = scaling(_df)
  _df = sign(_df)
  _df['d_et'] = _df['scaling']*_df['sign']
  return _df

def all_diagnostics(_df):
  """calcualtes all diagnostics"""
  _df = uwue(_df)
  _df = clean_df(_df)
  _df = pm_et(_df)
  _df = d_et_d_vpd(_df)
  return _df
