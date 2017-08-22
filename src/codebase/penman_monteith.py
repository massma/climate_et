#! ~/edward/bin/python
"""
this module has various functions for calculating
ET using penman-monteith
"""

import numpy as np
import pandas as pd
import metcalcs as met
from scipy.optimize import fsolve

#geophysical constats
VP_FACTOR = 100. #convert hPa -> Pa
K = 0.41 # vonkarmans constant
CP = 1012. # specific heat air, got from Changjie (I usually use 1004)
GAMMA = 66. # psychrometric constant
LV = 2.5e6
R_AIR = .0289645 # mean kg/mol dry air
R_STAR = 8.3144598 #J /mol/K

# From Changjie's oren model for stomatal resistance
# see "Survey and synthesis of intra- and interspecific variation
# in stomatal sensitivity to vapour pressure deficit" - Oren
# Gs = G1+G2*ln(VPD) in mm/s
OREN = pd.read_csv('../dat/orens_model.csv')
# convert to m/s
OREN.iloc[:, 2:6] = OREN.iloc[:, 2:6]/1000.
OREN.index = OREN.PFT

#see comments in files below for info on source, etc.
WUE_MEDLYN = pd.read_csv('../dat/franks_et_al_table2.csv',\
                         comment='#', delimiter=',')
WUE_MEDLYN.index = WUE_MEDLYN.PFT

WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
                  comment='#', delimiter=',')
WUE.index = WUE.PFT

LAI = pd.read_csv('../dat/bonan_et_al_table4_lai.csv',\
                  comment='#', delimiter=',')
LAI.index = LAI.PFT

DIM = pd.read_csv('../dat/bonan_et_al_table5_dimensions.csv',\
                  comment='#', delimiter=',')
DIM.index = DIM.PFT

ROUGH = pd.read_csv('../dat/wrf_roughness_length.csv',
          comment='#', delimiter=',')
ROUGH.index = ROUGH.PFT

MEDLYN = pd.read_csv('../dat/changjie_medlyn.csv',\
           comment='#', delimiter=',')
#convert to m/s
MEDLYN.iloc[:, 2:6] = MEDLYN.iloc[:, 2:6]/1000.
MEDLYN.index = MEDLYN.PFT


ADAM_MEDLYN = pd.read_csv('../dat/adam_mm_s_medlyn.csv',\
           comment='#', delimiter=',')
#now I take care of below conversion in adam funciton, b/c units
# for this are acutall mol/m2/s
# #convert to m/s - note only need to convert g0 b/c of functional form
# ADAM_MEDLYN.iloc[:, 1:3] = ADAM_MEDLYN.iloc[:, 1:3]/1000.
ADAM_MEDLYN.index = ADAM_MEDLYN.PFT

ADAM_MEDLYN_ET = pd.read_csv('../dat/adam_mm_s_W_m2_medlyn.csv',\
           comment='#', delimiter=',')
#now I take care of below conversion in adam funciton, b/c units
# for this are acutall mol/m2/s
# #convert to m/s - note only need to convert g0 b/c of functional form
# ADAM_MEDLYN.iloc[:, 1:3] = ADAM_MEDLYN.iloc[:, 1:3]/1000.
ADAM_MEDLYN_ET.index = ADAM_MEDLYN_ET.PFT


LEUNING = pd.read_csv('../dat/changjie_leuning.csv',\
           comment='#', delimiter=',')
#convert to m/s
LEUNING.iloc[:, 2:6] = LEUNING.iloc[:, 2:6]/1000.
LEUNING.index = LEUNING.PFT

FITTED_M = pd.read_csv('../dat/changjie_fitted_m.csv',\
           comment='#', delimiter=',')
#convert to m/s
FITTED_M.iloc[:, 2:6] = FITTED_M.iloc[:, 2:6]/1000.
FITTED_M.index = FITTED_M.PFT


def oren_r_e(vpd, pft):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in s/m
  """
  g_1 = OREN.loc[pft].G1mean
  g_2 = OREN.loc[pft].G2mean
  return 1./(g_1 + g_2*np.log(vpd/1000.))

def medlyn_g_w(vpd, co2, rho, pft, _et):
  """
  returns leaf stomatal conductance in m/s
  this function, unlike the others, uses LAI, but this also
  creates problems
  et : W/m2
  pft : planf functional type
  vpd : Pa
  co2 : ppm
  """
  #convert g C -> mu mol C
  wue = WUE.loc[pft, 'u_wue_yearly']*1.e6/12.011
  # note bellow assumes that atmos co2 is same as leaf, might be bad
  _g1 = WUE_MEDLYN.loc[pft, 'g1M'] # note this sqrt(kPa)
  # below is units mol air / m2 / s
  g_w = 1.6*(1. + _g1/np.sqrt(vpd/1000.))*wue*_et/LV/np.sqrt(vpd/100.)/co2
  g_w = g_w*R_AIR/rho
  return g_w

def medlyn_r_e(vpd, pft, _et):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in s/m
  """
  #convert g C -> mu mol C
  wue = WUE.loc[pft, 'u_wue_yearly']*1.e6/12.011
  g_0 = MEDLYN.loc[pft].G0mean
  g_1 = MEDLYN.loc[pft].G1mean
  return 1./(g_0 + g_1/np.sqrt(vpd/1000.)*wue*_et/LV/np.sqrt(vpd/100.))

def adam_medlyn_r_e(vpd, t_a, _canopy, _et):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in s/m
  """
  #convert g C -> mu mol C
  wue = WUE.loc[_canopy['pft'], 'u_wue_yearly']*1.e6/12.011
  if 'g0' in _canopy:
    g_0 = _canopy['g0']
    g_1 = _canopy['g1']
  else:
    g_0 = ADAM_MEDLYN.loc[pft].g0_mean/1000.
    g_1 = ADAM_MEDLYN.loc[pft].g1_mean
  _canopy.loc['r_s'] = 1./(g_0*wue*_et/LV/np.sqrt(vpd/100.)\
                           *(1. + g_1/np.sqrt(vpd/1000.)))
  return _canopy

def et_adam_medlyn_r_e(vpd, _canopy):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in W/m2 * s/m
  """
  #convert g C -> mu mol C
  wue = WUE.loc[_canopy['pft'].iloc[0], 'u_wue_yearly']*1.e6/12.011
  if 'g0' in _canopy:
    g_0 = _canopy['g0']
    g_1 = _canopy['g1']
  else:
    g_0 = ADAM_MEDLYN_ET.loc[pft].g0_mean/1000.
    g_1 = ADAM_MEDLYN_ET.loc[pft].g1_mean
  print(wue)
  _canopy['r_s'] = 1./(g_0*wue/LV/np.sqrt(vpd/100.)\
                       *(1. + g_1/np.sqrt(vpd/1000.)))
  return _canopy


def leuning_r_e(vpd, pft, _et):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in s/m
  """
  #convert g C -> mu mol C
  wue = WUE.loc[pft, 'u_wue_yearly']*1.e6/12.011
  g_0 = LEUNING.loc[pft].G0mean
  g_1 = LEUNING.loc[pft].G1mean
  return 1./(g_0 + g_1/(vpd/1000.)*wue*_et/LV/np.sqrt(vpd/100.))

def fitted_m_r_e(vpd, pft, _et):
  """
  calculates ecosystem canopy resistance given vpd and plant functional type
  vpd : vapor pressure deficit in Pa
  pft : three letter plant functional type
  returns canopy resistance in s/m
  """
  #convert g C -> mu mol C
  wue = WUE.loc[pft, 'u_wue_yearly']*1.e6/12.011
  g_0 = FITTED_M.loc[pft].G0mean
  g_1 = FITTED_M.loc[pft].G1mean
  _m = FITTED_M.loc[pft].m_mean
  return 1./(g_0 + g_1/(vpd/1000.)**_m*wue*_et/LV/np.sqrt(vpd/100.))

def psim(ksi):
  """
  From Changjie, who adapted from Paulson 1970 (unsstablel ksi <0)
  and Beljaars and Holtstag 1991 (stable ksi >0)
  """
  unstable_idx = (ksi < 0.)
  _psim = np.ones(ksi.shape)*np.nan
  x =(1. - 16.*ksi[unstable_idx])**0.25
  a = (1. + x)/2.
  b = (1. + x**2)/2.
  _psim[unstable_idx] = np.log(a**2*b)-2.*np.arctan(x)+np.arcsin(1.)
  _ksi = ksi[~unstable_idx]
  _psim[~unstable_idx]  = -(0.7*_ksi+0.75*(_ksi-5./0.35)\
                            *np.exp(-0.35*_ksi)+3.75/0.35)
  return _psim

def psih(ksi):
  """
  From Changjie, who adapted from Paulson 1970 (unsstablel ksi <0)
  and Beljaars and Holtstag 1991 (stable ksi >0), note comment
  that stable could be broken?
  """
  unstable_idx = (ksi < 0.)
  _psih = np.ones(ksi.shape)*np.nan
  x =(1. - 16.*ksi[unstable_idx])**0.25
  a = (1. + x**2)/2.
  _psih[unstable_idx] = 2.*np.log(a)
  _ksi = ksi[~unstable_idx]
  print(_ksi[~np.isnan(_ksi)].size)
  _psih[~unstable_idx]  = -((1.+(2.*_ksi)/3.)**1.5+0.667*(_ksi-5./0.35)\
                            *np.exp(-(0.35*_ksi))+(0.667*5.)/0.35-1.)
  return _psih

def r_a(_atmos, _canopy):
  """
  returns atmospheric resistance in s/m,
  see pg 298, eq. 20.36  in Shuttleworth
  implicit assumption that z0=z0h=z0m
  """
  return np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0'])**2\
    /K**2/_atmos['u_z']


def corrected_r_a(_atmos, _canopy):
  """
  returns atmospehric resistsance in s/m, but requires vars ustar and
  heat flux (h) in _atmos, which might not be actually available many
  of times. only gets called by penman_monteith if ustar is in _atmos
  """
  _atmos['L'] = -_atmos['ustar']**3*CP*_atmos['rho_a']*(_atmos['t_a']+273.15)\
                /(0.41*9.7*_atmos['h'])
  _r_a = np.ones(_atmos['L'].shape)*np.nan
  neutral_idx = (_atmos['L'] == 0.)
  if _r_a[neutral_idx].size > 0.:
    print(_r_a[neutral_idx].size)
    _r_a[neutral_idx] = r_a(_atmos, _canopy)[_atmos['L'] == 0.]
  #for neutral no correction required
  # dimensionless stability param
  _atmos['ksi'] = (_canopy['zmeas']-_canopy['d'])/_atmos['L']
  _atmos['ksi'][neutral_idx] = np.nan
  _atmos['psim'] = psim(_atmos['ksi'])
  _atmos['psih'] = psih(_atmos['ksi'])
  #use log profiles
  _ksi0 = (_canopy['zmeas']-_canopy['d'])/_canopy['z0']
  _log_ksi0 = np.log(_ksi0)
  _ = (_log_ksi0-_atmos['psih'])*(_log_ksi0-_atmos['psim'])\
                       /(K**2*_atmos['u_z'])
  _r_a[~neutral_idx] = _[~neutral_idx]
  return _r_a

def delta(_atmos):
  """calculates delta in Pa/C"""
  return 4098.*_atmos['e_s']/((237.3+_atmos['t_a'])**2)

def gamma(_atmos):
  """calculates gamma in Pa/C from pressure and T"""
  return 1.e-6*CP*_atmos['p_a']/(0.622*(2.501-0.00236*_atmos['t_a']))

def rho_air(_atmos):
  """returns rho in kg/m3"""
  return (_atmos['p_a']/(287.05*(_atmos['t_a']+273.15)))\
    *(1.-(0.378*_atmos['e_s']/_atmos['p_a']*_atmos['rh']/100.))

def penman_monteith_prep(_atmos, _canopy):
  """
  does calculations on data sructures
  atmos and canopy requred by penman moneith
  returns updated data structures
  """
  #derived constants
  _atmos['e_s'] = met.vapor_pres(_atmos['t_a'])*VP_FACTOR
  _atmos['delta'] = delta(_atmos)
  _atmos['gamma'] = gamma(_atmos)
  _atmos['rho_a'] = rho_air(_atmos)

  # below allows us to provide vpd
  if 'vpd' not in _atmos:
    _atmos['vpd'] = _atmos['e_s']*(1.-_atmos['rh']/100.)

  if 'g_flux' not in _canopy:
    _canopy['g_flux'] = 0.05*_atmos['r_n'] # soil heat flux (W/m2)

  if 'z0' not in _canopy:
    _canopy['d'] = 2./3.*_canopy['height'] #zero plain displacement
    _canopy['z0'] = 0.1*_canopy['height'] # height moisture source/sink

  if 'zmeas' not in _canopy:
    _canopy['zmeas'] = 2.+_canopy['height'] # measurement height

  if 'r_a' not in _atmos:
    if 'ustar' in _atmos:
      _atmos['r_a'] = corrected_r_a(_atmos, _canopy)
    else:
      _atmos['r_a'] = r_a(_atmos, _canopy)

  # if 'r_s' not in _canopy:
  #   _canopy['r_s'] = oren_r_e(_atmos['vpd'], _canopy['pft'])

  if 'height' not in _canopy:
    try:
      _dim = DIM.loc[_canopy['pft']]
    except KeyError:
      # dimesnsions are roughly either forest or not, so in case
      # we don't have dim. data for a PFT use that distiction
      if _canopy['pft'][-1] == 'F':
        _dim = DIM.loc['DBF']
      else:
        _dim = DIM.loc['SH']
    _rough = ROUGH.loc[_canopy['pft'], 'z0']
    _canopy['z0'] = _rough
    _canopy['height'] = _rough/_dim.rough_len_factor
    _canopy['d'] = _canopy['height']*_dim.displacement_height
    _canopy['zmeas'] = 2.+_canopy['height'] # measurement height

  return _atmos, _canopy
  
def penman_monteith(_atmos, _canopy):
  """
  returns ET in W/m2
  _atmos :: dict of atmospheric vars
  _canopy :: class of _canopy vars
  """
  _atmos, _canopy = penman_monteith_prep(_atmos, _canopy)

  _et = (_atmos['delta']*\
       (_atmos['r_n']-_canopy['g_flux'])+\
         _atmos['rho_a']*CP*_atmos['vpd']/_atmos['r_a'])\
       /(_atmos['delta']+_atmos['gamma']*(1. + _canopy['r_s']/_atmos['r_a']))
  return _et

def penman_monteith_uwue(_atmos, _canopy):
  """
  returns ET in W/m2, different from penman_monteith in that uWUE effeciency
  is used to change the functional form of the problem
  _atmos :: dict of atmospheric vars
  _canopy :: class of _canopy vars
  """
  if 'vpd_leaf' in _atmos:
    vpd = _atmos['vpd_leaf']
  else:
    vpd = _atmos['vpd']

  _canopy = et_adam_medlyn_r_e(vpd, _canopy)
  _atmos, _canopy = penman_monteith_prep(_atmos, _canopy)
  _atmos.loc[_atmos['r_a'] < 0.1, 'r_a'] = np.nan
  _et = (_atmos['delta']*\
       (_atmos['r_n']-_canopy['g_flux'])+\
         (_atmos['rho_a']*CP*_atmos['vpd'] - _canopy['r_s'])/_atmos['r_a'])\
         /(_atmos['delta']+_atmos['gamma'])
  return _et

def optimizer_wrapper(_et, *env_vars):
  """
  solves for ET using uWUE, designed to be called by
  scipy.optmize.fsolve
  """
  _atmos, _canopy = env_vars
  if 'vpd_leaf' in _atmos:
    vpd = _atmos['vpd_leaf']
  else:
    vpd = _atmos['vpd']

  if _canopy['stomatal_model'] == 'medlyn':
    _canopy['r_s'] = medlyn_r_e(vpd, _canopy['pft'], _et)
  elif _canopy['stomatal_model'] == 'leuning':
    _canopy['r_s'] = leuning_r_e(vpd, _canopy['pft'], _et)
  elif _canopy['stomatal_model'] == 'fitted_m':
    _canopy['r_s'] = fitted_m_r_e(vpd, _canopy['pft'], _et)
  elif _canopy['stomatal_model'] == 'adam_medlyn':
    _canopy = adam_medlyn_r_e(vpd, _atmos['t_a'],\
                              _canopy, _et)
  elif _canopy['stomatal_model'] == 'medlyn_lai':
    _canopy['r_s'] = 1./(_canopy['lai']\
              *medlyn_g_w(vpd, _atmos['co2'], _atmos['rho_a'],\
                    _canopy['pft'], _et))
  else:
    print("ERROR!!! Neither r_s nor stomatal_model defined!")

  # if _et == 0.:
  #   f_out = 0.
  # else:
  #   f_out = penman_monteith(_atmos, _canopy) - _et
  f_out = penman_monteith(_atmos, _canopy) - _et
  return f_out

def recursive_penman_monteith(_atmos, _canopy, et0=1000., name='et'):
  """
  This module solves for ET using scipy optmize fsolve with WUE.
  This is a relatively simple function so should always converge.
  Optional argument et0 is the first guess for et to pass to solver.
  """

  result = pd.Series(data=np.ones(_atmos.shape[0])*np.nan,\
                     index=_atmos.index, name=name)

  for index in result.index:
    result.loc[index] = fsolve(optimizer_wrapper, et0,\
                           args=(_atmos.loc[index, :].copy(),\
                                 _canopy.loc[index, :].copy()))

  return result
