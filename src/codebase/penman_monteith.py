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

#below taken from Franks et al (with Berry)
MEDLYN = pd.read_csv('../dat/franks_et_al_table2.csv',\
                     comment='#', delimiter=',')
MEDLYN.index = MEDLYN.PFT

WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
                     comment='#', delimiter=',')
WUE.index = WUE.PFT


def medlyn_g_w(vpd, co2, rho, pft, _et):
    """
    returns leaf stomatal conductance in m/s
    et : W/m2
    pft : planf functional type
    vpd : Pa
    co2 : ppm
    """
    #convert g C -> mu mol C
    wue = WUE.loc[pft, 'u_wue_yearly']*1.e6/12.011
    # note bellow assumes that atmos co2 is same as leaf, might be bad
    _g1 = MEDLYN.loc[pft, 'g1M'] # note this sqrt(kPa)
    # below is units mol air / m2 / s
    g_w = 1.6*(1. + _g1/np.sqrt(vpd/1000.))*wue*_et/LV/np.sqrt(vpd/100.)/co2
    g_w = g_w*R_AIR/rho
    return g_w


def oren_r_l(vpd, pft):
    """
    calculates canopy resistance given vpd and plant functional type
    vpd : vapor pressure deficit in Pa
    pft : three letter plant functional type
    returns canopy resistance in s/m
    """
    g_1 = OREN.loc[pft].G1mean
    g_2 = OREN.loc[pft].G2mean
    return 1./(g_1 + g_2*np.log(vpd/1000.))

def r_a(_atmos, _canopy):
    """
    returns atmospheric resistance in s/m,
    see pg 298, eq. 20.36  in Shuttleworth
    """
    return (np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0m'])\
            *np.log((_canopy['zmeas']-_canopy['d'])/_canopy['z0h']))\
            /K**2/_atmos['u_z']

def r_s(_atmos, _canopy):
    """returns _canopy mean stomatal resistance in s/m"""
    return oren_r_l(_atmos['vpd'], _canopy['pft'])/_canopy['lai']

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

    _canopy['d'] = 2./3.*_canopy['height'] #zero plain displacement
    _canopy['z0m'] = 0.1*_canopy['height'] # height moisture source/sink
    _canopy['z0h'] = _canopy['z0m'] # height of heat source/sink
    _canopy['zmeas'] = 2.+_canopy['height'] # measurement height

    _r_a = r_a(_atmos, _canopy)
    if 'r_s' not in _atmos:
        _atmos['r_s'] = r_s(_atmos, _canopy)

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
    _atmos['r_s'] = 1./(_canopy['lai']\
                   *medlyn_g_w(_atmos['vpd'], _atmos['co2'], _atmos['rho_a'],\
                               _canopy['pft'], _et))
    return penman_monteith(_atmos, _canopy) - _et

def medlyn_penman_monteith(_atmos, _canopy, et0=300.):
    """
    This module solves for ET using scipy optmize fsolve with WUE.
    This is a relatively simple function so should always converge.
    Optional argument et0 is the first guess for et to pass to solver.
    """
    #et0 = np.ones(atmos['vpd'].shape)*et0
    _et = fsolve(optimizer_wrapper, et0, args=(_atmos, _canopy))
    return _et
