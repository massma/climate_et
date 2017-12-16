import metcalcs as met
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

def delta(_df):
  """calculates delta in Pa/C, from shuttleworth equation 2.18"""
  return 4098.*_df['e_s']/((237.3+_df['t_a'])**2)

def gamma(_df):
  """calculates gamma in Pa/C from pressure and T"""
  return 1.e-6*CP*_df['p_a']/(0.622*(2.501-0.00236*_df['t_a']))

def r_moist(_df):
  """calculates the gas constant of moist air"""
  return R_DRY/(1.-(0.378*_df['e_s']\
                    /_df['p_a']*_df['rh']/100.))

def rho_air(_df):
  """returns rho in kg/m3"""
  return _df['p_a']/(_df['r_moist']*(_df['t_a']+273.15))

def saturation_vapor_pressure(_df):
  """returns asturation vapor pressure"""
  return met.vapor_pres(_df['t_a'])*VP_FACTOR

def psim(ksi):
  """
  From Changjie, who adapted from Paulson 1970 (unsstablel ksi <0)
  and Beljaars and Holtstag 1991 (stable ksi >0)
  """
  unstable_idx = (ksi < 0.)
  _psim = np.ones(ksi.shape)*np.nan
  chik2 = np.sqrt(1. - 16.*ksi[unstable_idx])
  chik  = np.sqrt(chik2)
  _psim[unstable_idx] = 2.*np.log((1.+chik)*0.5) + np.log((1.+chik2)*0.5)\
                        -2.*np.arctan(chik)+0.5*np.pi

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
  chik2 = np.sqrt(1. - 16.*ksi[unstable_idx])
  _psih[unstable_idx] = 2.*np.log((1.+chik2)*0.5)
  _ksi = ksi[~unstable_idx]
  print(_ksi[~np.isnan(_ksi)].size)
  _psih[~unstable_idx]  = -((1.+(2.*_ksi)/3.)**1.5+0.667*(_ksi-5./0.35)\
                            *np.exp(-(0.35*_ksi))+(0.667*5.)/0.35-1.)
  return _psih

def r_a(_df):
  """
  returns atmospheric resistance in s/m,
  see pg 298, eq. 20.36  in Shuttleworth
  implicit assumption that z0=z0h=z0m - changjie says maybe reduce factor 10
  """
  return np.log((_df['zmeas']-_df['d'])/_df['z0'])**2\
    /(K**2*_df['u_z'])

def corrected_r_a(_df):
  """
  returns atmospehric resistsance in s/m, but requires vars ustar and
  heat flux (h) in _df, which might not be actually available many
  of times. should only gets called when ustar is in _df
  """
  ksit = 0.465 # point of continuity in the stability profiles for heat

  # checked below from shuttleworth pg 292, should be good except t_a
  # should be  virtual (as should flux).
  _df['L'] = -_df['ustar']**3*CP*_df['rho_a']*(_df['t_a']+273.15)\
                /(K*G*_df['sensible'])
  _r_a = np.ones(_df['L'].shape)*np.nan
  neutral_idx = (np.absolute(_df['L']) <= 1.e-4)
  if _r_a[neutral_idx].size > 0.:
    _r_a[neutral_idx] = _df['r_a_uncorrected'][neutral_idx]
  _df['ksi'] = (_df['zmeas']-_df['d'])/_df['L']
  _df['ksi'][neutral_idx] = np.nan
  #below could be way more efficient
  _ = K/(np.log(-ksit*_df['L']/_df['z0'])
                 - psih(-ksit*np.ones(_df['L'].shape))
                 + psih(_df['z0']/_df['L'])
                 + 0.8*((ksit)**(-0.333)-(-_df['ksi'])**(-0.333)))
  ksit_idx = (_df['ksi'] < -ksit)
  _r_a[ksit_idx & (~neutral_idx)] = _[ksit_idx & (~neutral_idx)]
  _   = K/(np.log((_df['zmeas']-_df['d'])/_df['z0'])
                - psih(_df['ksi'])
                + psih(_df['z0']/_df['L']))
  zero_idx = (_df['ksi'] < 0.)
  _r_a[zero_idx & (~ksit_idx) & (~neutral_idx)] = _[zero_idx & (~ksit_idx)\
                                                    & (~neutral_idx)]
  _   = K/(np.log((_df['zmeas']-_df['d'])/_df['z0'])\
           + 5.*_df['ksi']-5.*_df['z0']/_df['L'])
  one_idx = (_df['ksi'] <= 1.)
  _r_a[one_idx & (~zero_idx) & (~ksit_idx)\
       & (~neutral_idx)] = _[one_idx & (~zero_idx) & (~ksit_idx)\
                             & (~neutral_idx)]
  _   = K/(np.log(_df['L']/_df['z0']) + 5. - 5.*_df['z0']/_df['L']
                + (5.*np.log(_df['ksi'])+_df['ksi']-1.))
  _r_a[(~one_idx) & (~zero_idx) & (~ksit_idx)\
       & (~neutral_idx)] = _[(~one_idx) & (~zero_idx) & (~ksit_idx)\
                             & (~neutral_idx)]
  _r_a[_r_a < 1.e-3] = 1.e-3
  _r_a = 1./(_r_a*_df['ustar'])
  return _r_a

def generate_vars(_df):
  """
  does calculations on data sructures
  atmos and canopy requred by penman moneith
  returns updated data structures
  """
  #derived constants
  _df['e_s'] = saturation_vapor_pressure(_df)
  _df['delta'] = delta(_df)
  _df['gamma'] = gamma(_df)
  _df['r_moist'] = r_moist(_df)
  _df['rho_a'] = rho_air(_df)

  # below allows us to provide vpd
  _df['vpd'] = _df['e_s']*(1.-_df['rh']/100.)

  _df['d'] = 2./3.*_df['height'] #zero plain displacement
    # roughness height, changjie emailed saying he use 1/10 of this for m/h
  _df['z0'] = 0.1*_df['height']


  _df['r_a_uncorrected'] = r_a(_df)
  _df['r_a_corrected'] = np.nan
  _df.loc[~np.isnan(_df.ustar),\
          'r_a_corrected'] = corrected_r_a(_df.loc[~np.isnan(_df.ustar), :])
  _df['r_a'] = np.nanmax(np.vstack([_df['r_a_corrected'],\
                                           _df['r_a_uncorrected']]), axis=0)
  _df['g_a'] = 1./_df['r_a']
  _df['t_a_k'] = _df['t_a'] + 273.15
  if 'c_a' not in _df:
    raise ValueError('c_a not loaded!!!!')

  return _df
