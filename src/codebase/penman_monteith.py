#! ~/edward/bin/python
"""
this module has various functions for calculating
ET using penman-monteith
"""

import numpy as np
import pandas as pd
import metcalcs as met


#geophysical constats
VP_FACTOR = 100. #convert hPa -> Pa
K = 0.41 # vonkarmans constant
CP = 1004. # specific heat air
GAMMA = 66. # psychrometric constant

# From Changjie's oren model for stomatal resistance
# see "Survey and synthesis of intra- and interspecific variation
# in stomatal sensitivity to vapour pressure deficit" - Oren
# Gs = G1+G2*ln(VPD) in mm/s
OREN = pd.read_csv('../../dat/orens_model.csv')
# convert to m/s
OREN.iloc[:, 2:6] = OREN.iloc[:, 2:6]/1000.
OREN.index = OREN.PFTs

def r_l(vpd, pft):
    """
    calculates canopy resistance given vpd and plant functional type
    vpd : vapor pressure deficit in Pa
    pft : three letter plant functional type
    returns canopy resistance in s/m
    """
    g_1 = OREN.loc[pft].G1mean
    g_2 = OREN.loc[pft].G2mean
    return 1./(g_1 + g_2*np.log(vpd/1000.))

def r_a(atmos, canopy):
    """returns atmospheric resistance in s/m"""
    return (np.log((canopy['zmeas']-canopy['d'])/canopy['z0m'])\
            *np.log((canopy['zmeas']-canopy['d'])/canopy['z0h']))\
            /K**2/atmos['u_z']

def r_s(atmos, canopy):
    """returns stomatal resistance in s/m"""
    return r_l(atmos['vpd'], canopy['pft'])/canopy['lai']

def penman_monteith(atmos, canopy):
    """
    returns ET in W/m2
    atmos :: dict of atmospheric vars
    canopy :: class of canopy vars
    """
    #derived constants
    atmos['delta'] = (met.vapor_pres(atmos['t_a']+0.1)\
                     -met.vapor_pres(atmos['t_a']))/0.1*VP_FACTOR
    atmos['e_s'] = met.vapor_pres(atmos['t_a'])*VP_FACTOR
    atmos['vpd'] = atmos['e_s']*(1.-atmos['rh']/100.)
    atmos['g_flux'] = 0.05*atmos['r_n'] # soil heat flux (W/m2)

    canopy['d'] = 2./3.*canopy['height']
    canopy['z0m'] = 0.1*canopy['height']
    canopy['z0h'] = canopy['z0m']
    canopy['zmeas'] = 2.+canopy['height']

    _r_a = r_a(atmos, canopy)
    _r_s = r_s(atmos, canopy)

    _et = (atmos['delta']*\
           (atmos['r_n']-atmos['g_flux'])+atmos['rho_a']*CP*atmos['vpd']/_r_a)\
           /(atmos['delta']+(1.+GAMMA)*(_r_s/_r_a))
    return _et
