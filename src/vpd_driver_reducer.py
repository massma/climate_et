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
from codebase.penman_monteith import penman_monteith
import metcalcs as met

# penamn monteith global warming impact calculation
FONT = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

mpl.rc('font', **FONT)

plt.close('all')

def plot_results(result, atmos):
    """plot GW experiment results"""

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_wireframe(atmos['rh'], atmos['t_a'],\
    #                   result)
    ax = fig.add_subplot(111)
    # vmax = np.absolute([result.max(), result.min()]).max()
    vmax = 3.*result.std()
    color = ax.pcolormesh(atmos['rh'], atmos['t_a'], result, cmap='RdBu',\
                          vmin=-vmax, vmax=vmax)
    ax.set_xlabel('RH')
    ax.set_ylabel('T')
    cbar = plt.colorbar(color)
    cbar.set_label(r'$\frac{\partial ET}{\partial VPD}$'\
                   ' ($W m^{-2}$  $Pa^{-1}$)')
    # plt.savefig('%s/temp/vpd.png' % os.environ['PLOTS'])
    plt.show(block=False)
    return


def d_et_d_vpd(atmos, canopy, pert, var='vpd'):
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
    var :: variable key derivative is taken wrt
    """
    result = (penman_monteith(atmos, canopy)\
              -penman_monteith(pert, canopy))\
             /(atmos[var]-pert[var])
    plot_results(result, atmos)
    return result

if str(__name__) == "__main__":
    ATMOS = {}
    ATMOS['r_n'] = 300. #W/m2
    ATMOS['rho_a'] = 1.205 #density kg/m3
    ATMOS['t_a'] = np.linspace(15., 35.) # C
    ATMOS['rh'] = np.linspace(20., 99.) # rel humdidity

    ATMOS['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    ATMOS['vpd'] = met.vapor_pres(ATMOS['t_a'])*(1.-ATMOS['rh']/100.)*100.


    CANOPY = {}
    CANOPY['pft'] = 'EBF'
    CANOPY['height'] = 10. # plant heigh m
    CANOPY['lai'] = 1. # leaf area index pierre says 1 max feasible

    ATMOS['t_a'], ATMOS['rh'],  = np.meshgrid(ATMOS['t_a'], ATMOS['rh'])

    PERT = copy.deepcopy(ATMOS)
    PERT['vpd'] += 1. # add 1 hpa for perturbation



    RESULT = d_et_d_vpd(ATMOS, CANOPY, PERT)
