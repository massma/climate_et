

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

def calc_uwue_old(_df):
  """cals uwue given structure, with an explcit et_obs
  This was the onld structure before refactoring"""
  return - _df['gamma']*_df['c_a']\
                *np.sqrt(_df['vpd'])*_df['p_a']\
                /(_df['r_a']\
                  *(_df['et_obs']*(_df['gamma']+_df['delta'])\
                    -_df['delta']*(_df['r_n']-_df['g_flux'])\
                    -1./_df['r_a']*_df['rho_a']*CP*_df['vpd'])
                  *1.6*R_STAR*(273.15 + _df['t_a'])\
                  *(1.+_df['g1']/np.sqrt(_df['vpd'])))

def calc_uwue(_df):
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
