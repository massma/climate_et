#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import metcalcs as met


# penamn monteith global warming impact calculation
FONT = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **FONT)

plt.close('all')

#geophysical constats
VP_FACTOR = 100. #convert hPa -> Pa
K = 0.41 # vonkarmans constant
CP = 1004. # specific heat air
GAMMA = 66. # psychrometric constant

# From Changjie's oren model for stomatal resistance
# see "Survey and synthesis of intra- and interspecific variation
# in stomatal sensitivity to vapour pressure deficit" - Oren
# Gs = G1+G2*ln(VPD) in mm/s
OREN = pd.read_csv('../dat/orens_model.csv')
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
    return (np.log((canopy.zmeas-canopy.d)/canopy.z0m)\
            *np.log((canopy.zmeas-canopy.d)/canopy.z0h))/K**2/atmos.u_z

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
    _r_a = r_a(atmos, canopy)
    _r_s = r_s(atmos, canopy)

    _et = (atmos['delta']*\
           (atmos['r_n']-atmos['g_flux'])+atmos['rho_a']*CP*atmos['vpd']/_r_a)\
           /(atmos['delta']+(1.+GAMMA)*(_r_s/_r_a))
    return _et

def gw_experiment_wrapper(atmos, canopy, gw_pert):
    """
    wrapper for script function
    atmos :: class of atmospheric vars
    canopy :: class of canopy vars
    gw_pert :: dictionary of global warming perturbation
    """

    cntrl = penman_monteith(atmos, canopy)

    result = {}
    for i in range(1, len(gw_pert.keys())+1):
        for pair in itertools.combinations(gw_pert.keys(), i):
            exp_dict = atmos.copy()
            newkey = ''
            for key in pair:
                exp_dict[key] += gw_pert[key]
            newkey = '-'.join(pair)
            result[newkey] = penman_monteith(exp_dict, canopy) - cntrl

    # fig = plt.figure()
    # ax1 = fig.add_subplot(211)
    # ax1.plot(co2, delta_et_r_n, '-b', label=r'$\delta ET_{rn}$')
    # ax1.plot(co2, delta_et_warming, '-r', label=r'$\delta ET_{warming}$')
    # ax1.plot(co2, delta_et_warming_nostoma, '-g', \
    #          label=r'$\delta ET_{warming no stoma}$')
    # ax1.plot(co2, delta_et_warming_vpdonly, '-c', \
    #          label=r'$\delta ET_{vpd only}$')
    # ax1.plot(co2, delta_et_full, '-r', linewidth=2,\
    #          label=r'$\delta ET_{full}$')
    # plt.legend(loc='best', fontsize=8)
    # #set(gca,'FontSize',16)
    # ax1.set_xlabel('[CO_2] (ppm)')
    # ax1.tick_params(right=True)#,axis=ax1.yaxis)
    # ax1.set_ylabel(r'$\delta$ ET (W/m**2)')
    # ax2 = fig.add_subplot(212)
    # ax2.plot(co2, 100*(delta_globalwarming/delta-1), \
    #          '-r', label=r'$\delta$ (#)')
    # ax2.plot(co2, 100*(vpd_globalwarming/delta_globalwarming*(delta/vpd)-1),\
    #          '-b', label=r'VPD / $\delta$ (#)')
    # #set(gca,'FontSize',16)
    # ax2.set_xlabel('[CO_2] (ppm)')
    # ax2.tick_params(right=True)#,axis=ax2.yaxis)
    # plt.legend(loc='best', fontsize=8)
    # # ylim([0 6])
    # plt.savefig('%s/temp/penman.png' % os.environ['PLOTS'])
    # plt.show(block=False)
    return result

if str(__name__) == "__main__":
    ATMOSPHEREENV = {}
    ATMOSPHEREENV['r_n'] = 300. #W/m2
    ATMOSPHEREENV['rho_a'] = 1.205 #density kg/m3
    ATMOSPHEREENV['t_a'] = 25. # C
    ATMOSPHEREENV['rh'] = 70. # rel humdidity
    ATMOSPHEREENV['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    ATMOSPHEREENV['delta'] = (met.vapor_pres(ATMOSPHEREENV['t_a']+0.1)\
                     -met.vapor_pres(ATMOSPHEREENV['t_a']))/0.1*VP_FACTOR
    ATMOSPHEREENV['e_s'] = met.vapor_pres(ATMOSPHEREENV['t_a'])*VP_FACTOR
    ATMOSPHEREENV['vpd'] = ATMOSPHEREENV['e_s']*(1.-ATMOSPHEREENV['rh']/100.)
    ATMOSPHEREENV['g_flux'] = 0.05*ATMOSPHEREENV['r_n'] # soil heat flux (W/m2)

    CANOPYENV = {}
    CANOPYENV['pft'] = 'EBF'
    CANOPYENV['height'] = 10. # plant heigh m
    CANOPYENV['lai'] = 1. # leaf area index pierre says 1 max feasible
    CANOPYENV['d'] = 2./3.*CANOPYENV['height']
    CANOPYENV['z0m'] = 0.1*CANOPYENV['height']
    CANOPYENV['z0h'] = CANOPYENV['z0m']
    CANOPYENV['zmeas'] = 2.+CANOPYENV['height']

    GW_PERT = {}
    GW_PERT['t_a'] = np.linspace(0., 3.7)
    GW_PERT['r_n'] = np.linspace(0, 8.5)
    CO2 = np.linspace(350., 1200.)

    RESULT = gw_experiment_wrapper(ATMOSPHEREENV, CANOPYENV, GW_PERT)
    print('done')
