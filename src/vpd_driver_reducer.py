#! ~/edward/bin/python
"""
this script explores the whether VPD is a driver or reducer
of ET, addressing task 1 in the readme.
"""
import os
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
    returns canopy resistance in m/s
    """
    g_1 = OREN.loc[pft].G1mean
    g_2 = OREN.loc[pft].G2mean
    return 1./(g_1 + g_2*np.log(vpd/1000.))

class AtmosphereEnv():
    """
    Atmospheric variables:
    r_n :: net radiation (W/m2)
    rho_a :: mean air density
    t_a :: t air in C
    rh :: percent relative humidity
    u_z :: wind speed at measurement height (m/s)
    """
    # pylint: disable=too-many-instance-attributes
    # 9 is reasonable in this case.
    # pylint : disable=too-few-public-methods
    # probably going to add functions to these classes
    def __init__(self, r_n=300., rho_a=1.205, t_a=25., rh=70., u_z=2.):
        self.r_n = r_n
        self.rho_a = rho_a
        self.t_a = t_a
        self.rh = rh
        self.u_z = u_z
        self.delta = (met.vapor_pres(self.t_a+0.1)\
                 -met.vapor_pres(self.t_a))/0.1*VP_FACTOR
        self.e_s = met.vapor_pres(self.t_a)*VP_FACTOR
        self.vpd = self.e_s*(1.-self.rh/100.)
        self.g_flux = 0.05*self.r_n # soil heat flux (W/m2)


class CanopyEnv():
    """
    canopy variables:
    pft :: plant functional type
    height :: crop height (m)
    lai :: leaf area idex (pierre says 1 is max feasible)
    """
    # pylint : disable=too-few-public-methods
    # probably going to add functions to these classes

    def __init__(self, pft='EBF', height=10., lai=1.):
        self.pft = pft
        self.height = height
        self.lai = lai
        self.d = 2./3.*self.height
        self.z0m = 0.1*self.height
        self.z0h = self.z0m
        self.zmeas = 2.+self.height



def wrapper(atmos=AtmosphereEnv(),\
            canopy=CanopyEnv()):
    """
    wrapper for script function
    atmos :: class of atmospheric vars
    canopy :: class of canopy vars
    """

    # gw stuff I iwll remove, just in now to make sure I don't
    # break script
    co2 = np.arange(350., 1200.) #[ppm]
    delta_co2 = co2[-1]-co2[0]
    dt_dco2 = 3.7/delta_co2
    dr_n_dco2 = 8.5/delta_co2

    #derived constants
    r_a = (np.log((canopy.zmeas-canopy.d)/canopy.z0m)\
           *np.log((canopy.zmeas-canopy.d)/canopy.z0h))/K**2/atmos.u_z
    r_s = r_l(atmos.vpd, canopy.pft)/canopy.lai

    lambda_et_ref = (atmos.delta*(atmos.r_n-g_flux)+atmos.rho_a*CP*atmos.vpd/r_a)/\
                    (atmos.delta+(1.+GAMMA)*(r_s/r_a))
    lambda_et_r_n = (atmos.delta*(atmos.r_n+dr_n_dco2*(co2-co2[1])-g_flux)\
                     +atmos.rho_a*CP*atmos.vpd/r_a)\
                     /(atmos.delta+(1+GAMMA)*(r_s/r_a))

    t_globalwarming = atmos.t_a + dt_dco2*(co2-co2[1])
    delta_globalwarming = (met.vapor_pres(t_globalwarming+0.1)\
                           -met.vapor_pres(t_globalwarming))/0.1*VP_FACTOR
    vpd_globalwarming = met.vapor_pres(t_globalwarming)\
                        *(1.-atmos.rh/100.)*VP_FACTOR
    r_s_globalwarming = r_l(vpd_globalwarming, canopy.pft)/canopy.lai

    lambda_et_warming = (delta_globalwarming*(atmos.r_n-g_flux)\
                         +atmos.rho_a*CP*vpd_globalwarming/r_a)\
                       /(delta_globalwarming+(1+GAMMA)*(r_s_globalwarming/r_a))
    lambda_et_warming_nostoma = (delta_globalwarming*(atmos.r_n-g_flux)\
                                +atmos.rho_a*CP*vpd_globalwarming/r_a)\
                                /(delta_globalwarming+(1+GAMMA)*(r_s/r_a))
    lambda_et_warming_vpdonly = (delta*(atmos.r_n-g_flux)\
                                 +atmos.rho_a*CP*vpd_globalwarming/r_a)\
                               /(delta+(1.+GAMMA)*(r_s_globalwarming/r_a))
    lambda_et_full = (delta_globalwarming\
                      *(atmos.r_n+dr_n_dco2*(co2-co2[1])-g_flux)\
                     +atmos.rho_a*CP*vpd_globalwarming/r_a)/\
                     (delta_globalwarming+(1.+GAMMA)*(r_s_globalwarming/r_a))

    delta_et_r_n = lambda_et_r_n - lambda_et_ref
    delta_et_warming = lambda_et_warming - lambda_et_ref
    delta_et_warming_nostoma = lambda_et_warming_nostoma - lambda_et_ref
    delta_et_warming_vpdonly = lambda_et_warming_vpdonly - lambda_et_ref
    delta_et_full = lambda_et_full - lambda_et_ref

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(co2, delta_et_r_n, '-b', label=r'$\delta ET_{rn}$')
    ax1.plot(co2, delta_et_warming, '-r', label=r'$\delta ET_{warming}$')
    ax1.plot(co2, delta_et_warming_nostoma, '-g', \
             label=r'$\delta ET_{warming no stoma}$')
    ax1.plot(co2, delta_et_warming_vpdonly, '-c', \
             label=r'$\delta ET_{vpd only}$')
    ax1.plot(co2, delta_et_full, '-r', linewidth=2,\
             label=r'$\delta ET_{full}$')
    plt.legend(loc='best', fontsize=8)
    #set(gca,'FontSize',16)
    ax1.set_xlabel('[CO_2] (ppm)')
    ax1.tick_params(right=True)#,axis=ax1.yaxis)
    ax1.set_ylabel(r'$\delta$ ET (W/m**2)')
    ax2 = fig.add_subplot(212)
    ax2.plot(co2, 100*(delta_globalwarming/delta-1), \
             '-r', label=r'$\delta$ (#)')
    ax2.plot(co2, 100*(vpd_globalwarming/delta_globalwarming*(delta/vpd)-1),\
             '-b', label=r'VPD / $\delta$ (#)')
    #set(gca,'FontSize',16)
    ax2.set_xlabel('[CO_2] (ppm)')
    ax2.tick_params(right=True)#,axis=ax2.yaxis)
    plt.legend(loc='best', fontsize=8)
    # ylim([0 6])
    plt.savefig('%s/temp/penman.png' % os.environ['PLOTS'])
    plt.show(block=False)
    return

if str(__name__) == "__main__":
    wrapper()
    print('done')
