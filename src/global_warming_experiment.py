#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import os
import copy
import itertools
import importlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import codebase.penman_monteith as pm
import metcalcs as met

importlib.reload(pm)

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
    plt.tight_layout()
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
    ctrl = pm.recursive_penman_monteith(atmos, canopy)

    # be careful below b/c if altering a buch of vars then the
    # combinations could get out of control
    for i in range(1, len(gw_pert.keys())+1):
        for pair in itertools.combinations(gw_pert.keys(), i):
            exp_dict = copy.deepcopy(atmos)
            newkey = ''
            for key in pair:
                exp_dict[key] += gw_pert[key]
            newkey = '-'.join(pair)
            exp_dict['vpd'] = met.vapor_pres(exp_dict['t_a'])\
                              *(1.-exp_dict['rh']/100.)*100.
            result[newkey] = pm.recursive_penman_monteith(exp_dict, canopy)-ctrl

    plot_results(result)
    return result

if str(__name__) == "__main__":
    ATMOS = {}
    ATMOS['r_n'] = 300. #W/m2
    ATMOS['rho_a'] = 1.205 #density kg/m3
    ATMOS['t_a'] = 20. # C
    ATMOS['rh'] = 80. # rel humdidity
    ATMOS['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    ATMOS['co2'] = 350. #ppm
    ATMOS['vpd'] = met.vapor_pres(ATMOS['t_a'])*(1.-ATMOS['rh']/100.)*100.


    CANOPY = {}
    CANOPY['pft'] = 'DBF'
    CANOPY['stomatal_model'] = 'medlyn'

    GW_PERT = {}
    GW_PERT['t_a'] = np.linspace(0., 3.7)
    GW_PERT['r_n'] = np.linspace(0, 8.5)
    GW_PERT['co2'] = np.linspace(350., 1200.) - 350.
    RESULT = gw_experiment_wrapper(ATMOS, CANOPY, GW_PERT)
