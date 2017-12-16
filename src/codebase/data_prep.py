import metcalcs as met

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

def r_a(_atmos, _canopy):
  """
  returns atmospheric resistance in s/m,
  see pg 298, eq. 20.36  in Shuttleworth
  implicit assumption that z0=z0h=z0m - changjie says maybe reduce factor 10
  """
  return np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0'])**2\
    /(K**2*_atmos['u_z'])

def corrected_r_a(_atmos, _canopy):
  """
  returns atmospehric resistsance in s/m, but requires vars ustar and
  heat flux (h) in _atmos, which might not be actually available many
  of times. should only gets called when ustar is in _atmos
  """
  ksit = 0.465 # point of continuity in the stability profiles for heat

  # checked below from shuttleworth pg 292, should be good except t_a
  # should be  virtual (as should flux).
  _atmos['L'] = -_atmos['ustar']**3*CP*_atmos['rho_a']*(_atmos['t_a']+273.15)\
                /(K*G*_atmos['h'])
  _r_a = np.ones(_atmos['L'].shape)*np.nan
  neutral_idx = (np.absolute(_atmos['L']) <= 1.e-4)
  if _r_a[neutral_idx].size > 0.:
    _r_a[neutral_idx] = _atmos['r_a_uncorrected'][neutral_idx]
  _atmos['ksi'] = (_canopy['zmeas']-_canopy['d'])/_atmos['L']
  _atmos['ksi'][neutral_idx] = np.nan
  #below could be way more efficient
  _ = K/(np.log(-ksit*_atmos['L']/_canopy['z0'])
                 - psih(-ksit*np.ones(_atmos['L'].shape))
                 + psih(_canopy['z0']/_atmos['L'])
                 + 0.8*((ksit)**(-0.333)-(-_atmos['ksi'])**(-0.333)))
  ksit_idx = (_atmos['ksi'] < -ksit)
  _r_a[ksit_idx & (~neutral_idx)] = _[ksit_idx & (~neutral_idx)]
  _   = K/(np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0'])
                - psih(_atmos['ksi'])
                + psih(_canopy['z0']/_atmos['L']))
  zero_idx = (_atmos['ksi'] < 0.)
  _r_a[zero_idx & (~ksit_idx) & (~neutral_idx)] = _[zero_idx & (~ksit_idx)\
                                                    & (~neutral_idx)]
  _   = K/(np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0'])\
           + 5.*_atmos['ksi']-5.*_canopy['z0']/_atmos['L'])
  one_idx = (_atmos['ksi'] <= 1.)
  _r_a[one_idx & (~zero_idx) & (~ksit_idx)\
       & (~neutral_idx)] = _[one_idx & (~zero_idx) & (~ksit_idx)\
                             & (~neutral_idx)]
  _   = K/(np.log(_atmos['L']/_canopy['z0']) + 5. - 5.*_canopy['z0']/_atmos['L']
                + (5.*np.log(_atmos['ksi'])+_atmos['ksi']-1.))
  _r_a[(~one_idx) & (~zero_idx) & (~ksit_idx)\
       & (~neutral_idx)] = _[(~one_idx) & (~zero_idx) & (~ksit_idx)\
                             & (~neutral_idx)]
  _r_a[_r_a < 1.e-3] = 1.e-3
  _r_a = 1./(_r_a*_atmos['ustar'])
  return _r_a

def generate_vars():
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
