#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import os
import copy
import importlib
from collections import OrderedDict
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import codebase.penman_monteith as pm
import metcalcs as met

importlib.reload(pm)

# penamn monteith global warming impact calculation
FONT = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

mpl.rc('font', **FONT)

plt.close('all')

def plot_results(result, atmos, canopy):
    """plot GW experiment results"""
    nplots = len(result)
    fig = plt.figure()
    fig.set_figheight(fig.get_figheight()*nplots)

    _ax = []
    for i, key in enumerate(result):
        # print(i)
        # print(key)
        _ax.append(fig.add_subplot(nplots, 1, i+1))
        vmax = result[key].mean() + 3.*result[key].std()
        vmin = -vmax
        if i > 0:
            cmap = 'viridis'
            vmax = None
            vmin = None
        else:
            cmap = 'RdBu'
        color = _ax[i].pcolormesh(atmos['rh'], atmos['t_a'], result[key],\
                                 cmap=cmap, vmin=vmin, vmax=vmax)
        _ax[i].set_xlabel('RH')
        _ax[i].set_ylabel('T')
        _ax[i].set_title('PFT: %s; %s VPD Changing'\
                         % (canopy['pft'], str(key)))
        cbar = plt.colorbar(color)
        cbar.set_label(r'$\frac{\partial ET}{\partial VPD_{%s}}$'\
                       ' ($W m^{-2}$  $Pa^{-1}$)' % key)
    plt.tight_layout()
    plt.savefig('%s/climate_et/%s_vpd.png'\
                % (os.environ['PLOTS'], canopy['pft']))
    plt.show(block=False)
    return


def d_et_d_vpd(atmos, canopy, pert, dvar=1.):
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
    pert :: dictionary of atmos variables with changed
    var :: variable key derivative is taken wrt
    thresh :: this is a p/m min/max threshold for the output
             sometimes the solver struggles, so this hopes to fix
    dvar :: this is what the difference in the variable is (defaul 1 hPa)
    """

    result = (pm.medlyn_penman_monteith(pert, canopy)\
                       -pm.medlyn_penman_monteith(atmos, canopy))\
                       /dvar
    return result

def main():
    pfts = ['CRO', 'DBF', 'C3G', 'C4G', 'ENF', 'SH']
    nvars = 100
    for pft in pfts:
        time_start = time.time()
        atmos = {}
        atmos['r_n'] = 300. #w/m2
        atmos['rho_a'] = 1.205 #density kg/m3
        atmos['t_a'] = np.linspace(15., 35., nvars) # c
        atmos['rh'] = np.linspace(40., 95., nvars) # rel humdidity
        atmos['co2'] = 400.
        # atmos['co2'] = np.array([400.,800.])
        # atmos['t_a'], atmos['rh'], atmos['co2'] = \
        #         np.meshgrid(atmos['t_a'], atmos['rh'], atmos['co2'])
        atmos['t_a'], atmos['rh'], = np.meshgrid(atmos['t_a'], atmos['rh'])

        atmos['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
        atmos['vpd'] = met.vapor_pres(atmos['t_a'])*(1.-atmos['rh']/100.)*100.
        #atmos['vpd_leaf'] = atmos['vpd'].copy()

        canopy = {}
        canopy['pft'] = pft
        canopy['lai'] = 1. # leaf area index pierre says 1 max feasible

        pert = copy.deepcopy(atmos)

        result = OrderedDict()

        pert['vpd_leaf'] = pert['vpd'] + 1.
        pert['vpd'] += 1. # add 1 pa for perturbation
        result['full'] = d_et_d_vpd(atmos, canopy, pert)

        pert['vpd'] -= 1.
        result['leaf'] = d_et_d_vpd(atmos, canopy, pert)

        pert['vpd_leaf'] -= 1.
        pert['vpd'] += 1.
        result['atm'] = d_et_d_vpd(atmos, canopy, pert)

        plot_results(result, atmos, canopy)
        print('for pft %s, time was %f s' % (pft, time.time()-time_start))
    return

if str(__name__) == "__main__":
    main()
    # TIME = time.time()
    # NVARS = 10
    # ATMOS = {}
    # ATMOS['r_n'] = 300. #W/m2
    # ATMOS['rho_a'] = 1.205 #density kg/m3
    # ATMOS['t_a'] = np.linspace(15., 35., NVARS) # C
    # ATMOS['rh'] = np.linspace(40., 95., NVARS) # rel humdidity
    # ATMOS['co2'] = 400.
    # # ATMOS['co2'] = np.array([400.,800.])
    # # ATMOS['t_a'], ATMOS['rh'], ATMOS['co2'] = \
    # #         np.meshgrid(ATMOS['t_a'], ATMOS['rh'], ATMOS['co2'])
    # ATMOS['t_a'], ATMOS['rh'], = np.meshgrid(ATMOS['t_a'], ATMOS['rh'])

    # ATMOS['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
    # ATMOS['vpd'] = met.vapor_pres(ATMOS['t_a'])*(1.-ATMOS['rh']/100.)*100.
    # #ATMOS['vpd_leaf'] = ATMOS['vpd'].copy()


    # CANOPY = {}
    # CANOPY['pft'] = 'DBF'
    # CANOPY['lai'] = 1. # leaf area index pierre says 1 max feasible


    # PERT = copy.deepcopy(ATMOS)

    # RESULT = OrderedDict()

    # PERT['vpd_leaf'] = PERT['vpd'] + 1.
    # PERT['vpd'] += 1. # add 1 PA for perturbation
    # RESULT['full'] = d_et_d_vpd(ATMOS, CANOPY, PERT)

    # PERT['vpd'] -= 1.
    # RESULT['leaf'] = d_et_d_vpd(ATMOS, CANOPY, PERT)

    # PERT['vpd_leaf'] -= 1.
    # PERT['vpd'] += 1.
    # RESULT['atm'] = d_et_d_vpd(ATMOS, CANOPY, PERT)

    # plot_results(RESULT, ATMOS)
    # print('time was %f s' % (time.time()-TIME))
