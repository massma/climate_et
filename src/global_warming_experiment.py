#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import os
import copy
import itertools
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from codebase.penman_monteith import penman_monteith

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


def gw_experiment_wrapper(atmos, canopy, gw_pert):
    """
    wrapper for doing simple GW experiments

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
    gw_pert :: dictionary of global warming perturbation
               vars and their change from 350 ppm conditions
    """

    result = {}
    ctrl = penman_monteith(atmos, canopy)

    # be careful below b/c if altering a buch of vars then the
    # combinations could get out of control
    for i in range(1, len(gw_pert.keys())+1):
        for pair in itertools.combinations(gw_pert.keys(), i):
            exp_dict = copy.deepcopy(atmos)
            newkey = ''
            for key in pair:
                exp_dict[key] += gw_pert[key]
            newkey = '-'.join(pair)
            result[newkey] = penman_monteith(exp_dict, canopy)-ctrl

    plot_results(result)
    return result

if str(__name__) == "__main__":
    ATMOSPHEREENV = {}
    ATMOSPHEREENV['r_n'] = 300. #W/m2
    ATMOSPHEREENV['rho_a'] = 1.205 #density kg/m3
    ATMOSPHEREENV['t_a'] = 25. # C
    ATMOSPHEREENV['rh'] = 70. # rel humdidity
    ATMOSPHEREENV['u_z'] = 2. #wind speed at meas. hiehgt (m/s)

    CANOPYENV = {}
    CANOPYENV['pft'] = 'EBF'
    CANOPYENV['height'] = 10. # plant heigh m
    CANOPYENV['lai'] = 1. # leaf area index pierre says 1 max feasible

    GW_PERT = {}
    GW_PERT['t_a'] = np.linspace(0., 3.7)
    GW_PERT['r_n'] = np.linspace(0, 8.5)
    RESULT = gw_experiment_wrapper(ATMOSPHEREENV, CANOPYENV, GW_PERT)

