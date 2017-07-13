#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import os
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from importlib import reload
import codebase.penman_monteith as penman_monteith
import metcalcs as met
reload(penman_monteith)

# penamn monteith global warming impact calculation
FONT = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **FONT)

plt.close('all')


def plot_results(result, co2max=1200.):
    """plot GW experiment results"""
    co2 = np.linspace(350., co2max)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for key in result:
        if co2.size != result[key].size:
            result[key] *= np.ones(co2.shape)
        ax1.plot(co2, result[key], label=r'$\delta$ %s' % (key))

    plt.legend(loc='best', fontsize=8)
    #set(gca,'FontSize',16)
    ax1.set_xlabel('[CO_2] (ppm)')
    ax1.tick_params(right=True)#,axis=ax1.yaxis)
    ax1.set_ylabel(r'$\delta$ ET (W/m**2)')
    plt.legend(loc='best', fontsize=8)
    plt.savefig('%s/temp/penman.png' % os.environ['PLOTS'])
    plt.show(block=False)
    return


def d_et_d_vpd(_atmos, canopy, pert):
    """
    numerically calculate dET/dpert

    atmos :: dict of atmospheric vars, should contain:
        r_n :: net radiation (W/m2)
        rho_a :: mean air density (kg/m3)
        t_a :: t air in C
        rh :: percent relative humidity
        u_z :: wind speed at measurement height (m/s)

    canopy :: dict of canopy vars, should contain:
        pft :: plant functional type
        height :: crop height (m)
        lai :: leaf area index (pierre says 1 is max feasible)
    pert :: dictionary of atmos variables with 1 changed
    """
    result = (penman_monteith.penman_monteith(_atmos, canopy)\
              -penman_monteith.penman_monteith(pert, canopy))\
             /(_atmos['vpd']-pert['vpd'])

    return result

if str(__name__) == "__main__":
    ATMOS = {}
    ATMOS['r_n'] = 300. #W/m2
    ATMOS['rho_a'] = 1.205 #density kg/m3
    ATMOS['t_a'] = t_a = np.linspace(0., 35.) # C
    ATMOS['rh'] = np.linspace(20., 99.) # rel humdidity
    ATMOS['t_a'], ATMOS['rh'] = np.meshgrid(ATMOS['t_a'], ATMOS['rh'])
    ATMOS['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    ATMOS['vpd'] = met.vapor_pres(ATMOS['t_a'])*(1.-ATMOS['rh']/100.)*100.

    print(ATMOS['vpd'].max())
    PERT = copy.deepcopy(ATMOS)
    PERT['vpd'] += 1. # add 1 hpa for perturbation
    print(ATMOS['vpd'].max())
    print(PERT['vpd'].max())
          
    
    CANOPY = {}
    CANOPY['pft'] = 'EBF'
    CANOPY['height'] = 10. # plant heigh m
    CANOPY['lai'] = 1. # leaf area index pierre says 1 max feasible

    RESULT = d_et_d_vpd(ATMOS, CANOPY, PERT)
