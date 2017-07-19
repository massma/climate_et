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
import codebase.penman_monteith as pm
import metcalcs as met

import importlib

importlib.reload(pm)

# penamn monteith global warming impact calculation
FONT = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

mpl.rc('font', **FONT)

plt.close('all')

def plot_results(result, atmos):
    """plot GW experiment results"""

    fig = plt.figure()
    fig.set_figheight(fig.get_figheight()*2)
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_wireframe(atmos['rh'], atmos['t_a'],\
    #                   result)
    ax1 = fig.add_subplot(211)
    # vmax = np.absolute([result.max(), result.min()]).max()
    vmax = 3.*result[:,:,0].std()
    color = ax1.pcolormesh(atmos['rh'][:,:,0], atmos['t_a'][:,:,0],\
                           result[:,:,0], cmap='RdBu', vmin=-vmax, vmax=vmax)
    ax1.set_xlabel('RH')
    ax1.set_ylabel('T')
    ax1.set_title('CO2: %f' % atmos['co2'][:,:,0].mean())
    cbar = plt.colorbar(color)
    cbar.set_label(r'$\frac{\partial ET}{\partial VPD}$'\
                   ' ($W m^{-2}$  $Pa^{-1}$)')

    vmax = 3.*result[:,:,1].std()
    ax2 = fig.add_subplot(212)
    color = ax2.pcolormesh(atmos['rh'][:,:,1], atmos['t_a'][:,:,1],\
                           result[:,:,1], cmap='RdBu', vmin=-vmax, vmax=vmax)
    ax2.set_xlabel('RH')
    ax2.set_ylabel('T')
    ax2.set_title('CO2: %f' % atmos['co2'][:,:,1].mean())
    cbar = plt.colorbar(color)
    cbar.set_label(r'$\frac{\partial ET}{\partial VPD}$'\
                   ' ($W m^{-2}$  $Pa^{-1}$)')

    plt.savefig('%s/temp/vpd.png' % os.environ['PLOTS'])
    plt.show(block=False)
    return


def d_et_d_vpd(atmos, canopy, pert, var='vpd', thresh=2.):
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
    thresh :: this is a p/m min/max threshold for the output
             sometimes the solver struggles, so this hopes to fix
    """
    # there is probably a much better way to do below
    _atmos = atmos.copy()
    _canopy = canopy.copy()
    it = np.nditer([atmos['t_a'], atmos['rh'], atmos['vpd'],\
                    atmos['co2'], pert[var], None])
    for t_a, rh, vpd, co2, pert_vpd, result in it:
        _atmos['t_a'] = t_a
        _atmos['rh'] = rh
        _atmos['vpd'] = vpd
        _atmos['co2'] = co2
        _pert = _atmos.copy()
        _pert[var] = pert_vpd
        result[...] = (pm.medlyn_penman_monteith(_atmos, _canopy)\
                       -pm.medlyn_penman_monteith(_pert, _canopy))\
                       /(_atmos[var]-_pert[var])


    results = it.operands[-1]
    # results[results < -thresh] = np.nan
    # results[results > thresh] = np.nan
    plot_results(results, atmos)
    return results

if str(__name__) == "__main__":
    ATMOS = {}
    ATMOS['r_n'] = 300. #W/m2
    ATMOS['rho_a'] = 1.205 #density kg/m3
    ATMOS['t_a'] = np.linspace(15., 35.,50) # C
    ATMOS['rh'] = np.linspace(40., 95.,50) # rel humdidity
    ATMOS['co2'] = np.array([400.,800.])

    ATMOS['t_a'], ATMOS['rh'], ATMOS['co2'] = \
            np.meshgrid(ATMOS['t_a'], ATMOS['rh'], ATMOS['co2'])

    ATMOS['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    ATMOS['vpd'] = met.vapor_pres(ATMOS['t_a'])*(1.-ATMOS['rh']/100.)*100.


    CANOPY = {}
    CANOPY['pft'] = 'DBF'
    CANOPY['height'] = 10. # plant heigh m
    CANOPY['lai'] = 1. # leaf area index pierre says 1 max feasible


    PERT = copy.deepcopy(ATMOS)
    PERT['vpd'] += 1. # add 1 PA for perturbation

    result = d_et_d_vpd(ATMOS, CANOPY, PERT)
