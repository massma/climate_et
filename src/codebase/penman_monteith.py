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
CP = 1004. # specific heat air
GAMMA = 66. # psychrometric constant
LV = 2.5e6
R_AIR = .0289645 # mean kg/mol dry air

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

def r_a(_atmos, _canopy):
    """
    returns atmospheric resistance in s/m,
    see pg 298, eq. 20.36  in Shuttleworth
    """
    return (np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0m'])\
            *np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0h']))\
            /K**2/_atmos['u_z']

def penman_monteith(_atmos, _canopy):
    """
    returns ET in W/m2
    _atmos :: dict of atmospheric vars
    _canopy :: class of _canopy vars
    """
    #derived constants
    _atmos['delta'] = (met.vapor_pres(_atmos['t_a']+0.1)\
                     -met.vapor_pres(_atmos['t_a']))/0.1*VP_FACTOR
    _atmos['e_s'] = met.vapor_pres(_atmos['t_a'])*VP_FACTOR

    # below allows us to provide vpd
    if 'vpd' not in _atmos:
        _atmos['vpd'] = _atmos['e_s']*(1.-_atmos['rh']/100.)
    _atmos['g_flux'] = 0.05*_atmos['r_n'] # soil heat flux (W/m2)

    if 'z0m' not in _canopy:
        _canopy['d'] = 2./3.*_canopy['height'] #zero plain displacement
        _canopy['z0m'] = 0.1*_canopy['height'] # height moisture source/sink
        _canopy['z0h'] = _canopy['z0m'] # height of heat source/sink
        _canopy['zmeas'] = 2.+_canopy['height'] # measurement height

    _r_a = r_a(_atmos, _canopy)
    if 'r_s' not in _atmos:
        _atmos['r_s'] = oren_r_e(_atmos['vpd'], _canopy['pft'])

    _et = (_atmos['delta']*\
           (_atmos['r_n']-_atmos['g_flux'])+\
           _atmos['rho_a']*CP*_atmos['vpd']/_r_a)\
           /(_atmos['delta']+GAMMA*(1. + _atmos['r_s']/_r_a))
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
        _atmos['r_s'] = medlyn_r_e(vpd, _canopy['pft'], _et)
    elif _canopy['stomatal_model'] == 'leuning':
        _atmos['r_s'] = leuning_r_e(vpd, _canopy['pft'], _et)
    elif _canopy['stomatal_model'] == 'fitted_m':
        _atmos['r_s'] = fitted_m_r_e(vpd, _canopy['pft'], _et)
    elif _canopy['stomatal_model'] == 'medlyn_lai':
        _atmos['r_s'] = 1./(_canopy['lai']\
                            *medlyn_g_w(vpd, _atmos['co2'], _atmos['rho_a'],\
                                        _canopy['pft'], _et))
    else:
        print("ERROR!!! Neither r_s nor stomatal_model defined!")

    # if _et == 0.:
    #     f_out = 0.
    # else:
    #     f_out = penman_monteith(_atmos, _canopy) - _et
    f_out = penman_monteith(_atmos, _canopy) - _et
    return f_out

def recursive_penman_monteith(_atmos, _canopy, et0=1000.):
    """
    This module solves for ET using scipy optmize fsolve with WUE.
    This is a relatively simple function so should always converge.
    Optional argument et0 is the first guess for et to pass to solver.
    """
    #et0 = np.ones(atmos['vpd'].shape)*et0
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
        _canopy['z0m'] = _rough
        _canopy['height'] = _rough/_dim.rough_len_factor
        _canopy['d'] = _canopy['height']*_dim.displacement_height
        _canopy['z0h'] = _canopy['z0m'] # height of heat source/sink
        _canopy['zmeas'] = 2.+_canopy['height'] # measurement height

    packed_array = []
    keys = _atmos.keys()
    for key in _atmos:
        packed_array.append(np.copy(_atmos[key]))
    packed_array.append(None)

    _it = np.nditer(packed_array)
    for array in _it:
        _atmos = {}
        for i, key in enumerate(keys):
            _atmos[key] = array[i]
        result = array[-1]
        result[...] = fsolve(optimizer_wrapper, et0, args=(_atmos, _canopy))

    return _it.operands[-1]
