#! ~/edward/bin/python

"""
This module calculates all diagnostics used in the paper
"""
import time
import numpy as np
import codebase.data_prep as d_prep

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

def medlyn_kernel(_df, vpd, g1=None):
  """calcs medlyn without g/c_s/1.6 multipliers"""
  if g1 is None:
    g1 = _df['g1']
  return 1.0+g1/np.sqrt(vpd)

def bb_kernel(_df, vpd, g1=None):
  """calcs kernel for ball-berry model"""
  if g1 is None:
    g1 = _df['g1b']
  return (g1*(1.0-vpd/_df['e_s']))/1.6

def uwue(_df, n=0.5, kernel=medlyn_kernel, g1=None):
  """
  calcs uwue from _df with obs
  note approx using c_a instead of c_s
  """
  return -_df['g_a']*_df['gamma']*_df['c_a']\
    *_df['vpd']**n*_df['p_a']\
    /((_df['et_obs']*(_df['delta']+_df['gamma'])\
       -_df['delta']*_df['r_net']\
       -_df['g_a']*_df['rho_a']*CP*_df['vpd'])\
      *1.6*R_STAR*_df['t_a_k']*(kernel(_df, vpd=_df['vpd'], g1=g1)))

def g1_medlyn(_df, uwue=None):
  """
  calcs g1 slope coefficent given data,
  note currently only works with medlyn model
  rather than all others
  """
  if uwue is None:
    uwue = _df['uwue']
  return (-_df['g_a']*_df['gamma']*_df['c_a']\
    *np.sqrt(_df['vpd'])*_df['p_a']\
    /((_df['et_obs']*(_df['delta']+_df['gamma'])\
       -_df['delta']*_df['r_net']\
       -_df['g_a']*_df['rho_a']*CP*_df['vpd'])\
      *1.6*R_STAR*_df['t_a_k']*uwue)\
      -1.0)*np.sqrt(_df['vpd'])

def solve(initial_uwue, initial_g1, _df, convergence=0.01, max_iter=1000):
  """claculates optimal uwue and g1, to some covergence,
  which is given as a fraction of the intiial guess"""
  print("for %s initial uwue guess is %f" % (_df.pft.iloc[0], initial_uwue))
  print("for %s initial g1 guess is %f" % (_df.pft.iloc[0], initial_g1))
  old_uwue_guess = initial_uwue
  old_g1_guess = initial_g1
  uwue_converge = convergence*initial_uwue
  g1_converge = convergence*initial_g1
  for i in range(max_iter):

    uwue_s = uwue(_df, g1=old_g1_guess)
    uwue_guess = np.nanmean(uwue_s)
    g1_s = g1_medlyn(_df, uwue=np.mean([old_uwue_guess, uwue_guess]))
    g1_guess = np.nanmean(g1_s)

    # print("for iter %d uwue guess is %f" % (i, uwue_guess))
    # print("for iter %d  g1 guess is %f" % (i, g1_guess))
    if ((np.absolute(g1_guess-old_g1_guess) < g1_converge) & \
        (np.absolute(uwue_guess-old_uwue_guess) < uwue_converge)):
      print("Success! uwue guess is %f" % (uwue_guess))
      print("Success! g1 guess is %f" % (g1_guess))
      return np.mean([uwue_guess, old_uwue_guess]),\
        np.mean([g1_guess, old_g1_guess])
    old_uwue_guess = np.mean([uwue_guess, old_uwue_guess])
    old_g1_guess = np.mean([g1_guess, old_g1_guess])
  #raise ValueError("Solution did not converge")
  print("error, solution did not converge, returning orig vals anyways")
  print("guesses were: uwue: %f, g1: %f" % (uwue_guess, g1_guess))
  return uwue_guess, g1_guess

def pft_solve_wrapper(_df):
  """designed to be called on groupby pft, calls solve"""
  print("\n\n\n **** starting to work on pft %s ****" % _df.pft.iloc[0])
  start = time.time()
  _df['g1_lin'] = _df['g1'].copy()
  uwue, g1 = solve(_df.uwue_zhou.iloc[0], _df.g1.iloc[0], _df)
  print("time for pft %s was %d s" % (_df.pft.iloc[0], (time.time()-start)))
  _df['uwue_fit'] = uwue
  _df['g1'] = g1
  return _df

def clean_df(_df, var='uwue', vpd_thresh=10.0):
  """
  cleans df, for now removing lower and upper 5 percent of values,
  based on calc'd uwue. Could look into more sophisticated outlier
  identification (seems like from histogram there are clear ouliers).
  """
  cleaned_df = _df.loc[((_df[var] < _df[var].quantile(q=0.95)) &\
                      (_df[var] > _df[var].quantile(q=0.05)) &\
                      (_df.vpd > vpd_thresh))].copy()
  return cleaned_df


def pm_et(_df, vpd=None, uwue=None, n=0.5, kernel=medlyn_kernel):
  """calculates et using our new penman-monteith"""
  if vpd is None:
    vpd = _df['vpd']
  if uwue is None:
    uwue = _df['uwue']
  return (_df['delta']*_df['r_net']\
    +_df['g_a']*_df['p_a']/_df['t_a_k']\
    *(CP*vpd/_df['r_moist']\
      -_df['gamma']*_df['c_a']*vpd**n\
      /(R_STAR*1.6*uwue\
        *(kernel(_df, vpd)))))\
        /(_df['delta']+_df['gamma'])

def sign(_df, vpd=None, uwue=None, g1=None, n=0.5):
  """
  calculates the 'sign' term of et derivative w.r.t. vpd,
  as a function of n (exponent in *WUE metric)
  """
  if vpd is None:
    vpd = _df['vpd']
  if uwue is None:
    uwue = _df['uwue']
  if g1 is None:
    g1 = _df['g1']
  return CP/_df['r_moist']\
                -_df['gamma']*_df['c_a']\
                /(1.6*R_STAR*uwue)\
                *(vpd**(n-1.0)*((n+0.5)*g1/np.sqrt(vpd)+n))\
                  /(g1/np.sqrt(vpd)+1.0)**2

def scaling(_df, t_a=None, g_a=None):
  """
  calculates the 'scaling' term of the et derivative w.r.t. vpd,
  note 2.0 factor is from esat-rh decomposition,
  but need to double check this
  """
  if t_a is None:
    delta = _df['delta']
  else:
    delta = d_prep.delta({'e_s' : d_prep.sat_vapor_press({'t_a' : t_a}),\
                          't_a' : t_a})
  if g_a is None:
    g_a = _df['g_a']
  return 2.0*g_a*_df['p_a']\
    /(_df['t_a_k']*(delta + _df['gamma']))


def medlyn(_df, vpd=None, gpp=None):
  """calcs medlyn given gpp obs"""
  if gpp is None:
    gpp = _df['gpp_obs']
  return R_STAR*_df['t_a_k']/_df['p_a']\
    *1.6*(1.0 + _df['g1']/np.sqrt(vpd))*gpp/_df['c_a']

def lai(_df, vpd=None, gpp=None):
  """calcs a lai given gpp obs"""
  if vpd is None:
    vpd = _df['vpd']
  if gpp is None:
    gpp = _df['gpp_obs']
  return _df['g_a']/(((_df['delta']*_df['r_net']\
                       +_df['g_a']*_df['p_a']*CP*vpd\
                       /(_df['t_a_k']*_df['r_moist']))\
                      /_df['et_obs']\
                      -_df['delta'])/_df['gamma'] -1.0)\
                      /medlyn(_df, vpd=_df['vpd'], gpp=gpp)

def pm_et_orig(_df, vpd=None, gpp=None, lai=None):
  """original penman monteith as a function of GPP"""
  if vpd is None:
    vpd = _df['vpd']
  if gpp is None:
    gpp = _df['gpp_obs']
  if lai is None:
    lai = _df['lai_pm']
  return (_df['delta']*_df['r_net']\
          +_df['g_a']*_df['p_a']*CP*vpd/(_df['t_a_k']*_df['r_moist']))\
          /(_df['delta']+_df['gamma']\
            *(1.0 + _df['g_a']/(lai*medlyn(_df, vpd=vpd, gpp=gpp))))

def et_leaf(_df, vpd=None, gpp=None, lai=None):
  """ecoystem as a leaf, no energy balance"""
  if vpd is None:
    vpd = _df['vpd']
  if gpp is None:
    gpp = _df['gpp_obs']
  if lai is None:
    lai = _df['lai_pm']
  return 1./(1./_df['g_a']+1./(lai*medlyn(_df, vpd=vpd, gpp=gpp)))\
    *_df['p_a']*CP*vpd/(_df['t_a_k']*_df['r_moist']*_df['gamma'])

def pet(_df, vpd=None):
  """caluclates pet"""
  if vpd is None:
    vpd = _df['vpd']
  return (_df['delta']*_df['r_net']\
          +_df['g_a']*_df['p_a']*CP*vpd/(_df['t_a_k']*_df['r_moist']))\
          /(_df['delta']+_df['gamma'])

def gpp_fixed_wrapper(_df, mean_df):
  """
  calculates et with original penman monteith, assuming
  GPP and lai is constant within PFT. used to compare with new
  penman monteith.
  Designed to be called by groupby on pft
  """
  _df['lai_gppfixed'] = lai(_df, gpp=mean_df.loc[_df.pft.iloc[0],\
                                                 'gpp_obs'])
  _df['et_gppfixed'] = pm_et_orig(_df, gpp=mean_df.loc[_df.pft.iloc[0],\
                                                       'gpp_obs'],\
                                  lai=_df['lai_gppfixed'].mean()) # 1.0
  return _df

def calc_sigma(_df):
  """calcualtes a sgima for variabiility"""
  # _df['sigma'] = _df['uwue']/_df['uwue_fit'].iloc[0]
  _df['sigma'] = _df['uwue']/_df['uwue'].mean()
  return _df

def all_diagnostics(_df):
  """calcualtes all diagnostics"""
  _df['uwue'] = uwue(_df)
  _df['uwue_bb'] = uwue(_df, kernel=bb_kernel)
  _df['iwue_bb'] = uwue(_df, kernel=bb_kernel, n=1.0)
  _df['iwue'] = uwue(_df, n=1.0)
  _df = _df.groupby('pft').apply(clean_df)
  _df = _df.reset_index(drop=True)
  # _df = _df.groupby('pft').apply(pft_solve_wrapper)
  # _df = _df.reset_index(drop=True)
  # recalc uwue with new fit g1
  # _df['uwue'] = uwue(_df)
  _df['et'] = pm_et(_df)
  _df['scaling'] = scaling(_df)
  _df['sign'] = sign(_df)
  _df['d_et'] = _df['scaling']*_df['sign']
  _df['lai_pm'] = lai(_df)
  _df['et_pm_original'] = pm_et_orig(_df)
  _df['pet'] = pet(_df)
  # no, do comparison to mean to set mean sigma = 1
  _df = _df.groupby('pft').apply(calc_sigma)
  _df = _df.reset_index(drop=True)
  dfs = {'mean' : _df.groupby('pft').mean(),\
         'min' : _df.groupby('pft').min(),\
         'median' : _df.groupby('pft').median(),\
         'max' : _df.groupby('pft').max(),\
         '5' : _df.groupby('pft').quantile(q=0.05),\
         '95' : _df.groupby('pft').quantile(q=0.95)}
  _df = _df.groupby('pft').apply(gpp_fixed_wrapper, dfs['mean'])
  _df = _df.reset_index(drop=True)
  dfs['full'] = _df
  return dfs
