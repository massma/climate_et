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
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **font)


plt.close('all')

# gw stuff I iwll remove, just in now to make sure I don't
# break script
CO2 = np.arange(350., 1200.) #[ppm]
DeltaCO2 = CO2[-1]-CO2[0]
dT_dCO2 = 3.7/DeltaCO2
dRn_dCO2 = 8.5/DeltaCO2

# From Changjie's oren model for stomatal resistance
# see "Survey and synthesis of intra- and interspecific variation
# in stomatal sensitivity to vapour pressure deficit" - Oren
# Gs = G1+G2*ln(VPD) in mm/s
oren = pd.read_csv('../dat/orens_model.csv')
# convert to m/s
oren.iloc[:, 2:6] = oren.iloc[:, 2:6]/1000.
oren.index = oren.PFTs

# user adjustable constants
PFT = 'EBF' #evergreen broadleaf
Rn = 300. # net radiation (W/m2)
G1 = oren.loc[PFT].G1mean
G2 = oren.loc[PFT].G2mean
rhoA = 1.205 # mean air density
Ta = 25. # C 273.15+
RH = 70. # precent
uz = 2. # wind speed at height
h = 10. # crop height
rlmin = 20.
LAI = 1. # max feasible LAI

#geophysical constats
vp_factor = 100. #convert hPa -> Pa
k = 0.41 # vonkarmans constant
cP = 1004. # specific heat air
gamma = 66. # psychrometric constant

#derived constants
G = 0.05*Rn # soil heat flux (W/m2)
d = 2./3.*h
z0m = 0.1*h
z0h = z0m
zmeas = 2.+h
ra = (np.log((zmeas-d)/z0m)*np.log((zmeas-d)/z0h))/k**2/uz
Delta = (met.vapor_pres(Ta+0.1)-met.vapor_pres(Ta))/0.1*vp_factor
es = met.vapor_pres(Ta)*vp_factor

VPD = es*(1.-RH/100.)
rl = 1./(G1 + G2*np.log(VPD/1000.)) #believe oren's take 
rs = rl/LAI

lambdaETref = (Delta*(Rn-G)+rhoA*cP*VPD/ra)/(Delta+(1.+gamma)*(rs/ra))
lambdaET_Rn = (Delta*(Rn+dRn_dCO2*(CO2-CO2[1])-G)+rhoA*cP*VPD/ra)\
              /(Delta+(1+gamma)*(rs/ra))

T_globalwarming = Ta + dT_dCO2*(CO2-CO2[1])
Delta_globalwarming = (met.vapor_pres(T_globalwarming+0.1)\
                       -met.vapor_pres(T_globalwarming))/0.1*vp_factor
VPD_globalwarming = met.vapor_pres(T_globalwarming)*(1.-RH/100.)*vp_factor
rs_globalwarming = 1/(G1+G2*np.log(VPD_globalwarming/1000.))/LAI

lambdaET_warming = (Delta_globalwarming*(Rn-G)+rhoA*cP*VPD_globalwarming/ra)\
                   /(Delta_globalwarming+(1+gamma)*(rs_globalwarming/ra))
lambdaET_warming_nostoma = (Delta_globalwarming*(Rn-G)\
                            +rhoA*cP*VPD_globalwarming/ra)\
                            /(Delta_globalwarming+(1+gamma)*(rs/ra))
lambdaET_warming_VPDonly = (Delta*(Rn-G)+rhoA*cP*VPD_globalwarming/ra)\
                           /(Delta+(1+gamma)*(rs_globalwarming/ra))
lambdaET_full = (Delta_globalwarming*(Rn+dRn_dCO2*(CO2-CO2[1])-G)\
                 +rhoA*cP*VPD_globalwarming/ra)/\
                 (Delta_globalwarming+(1+gamma)*(rs_globalwarming/ra))

DeltaET_Rn = lambdaET_Rn - lambdaETref
DeltaET_warming = lambdaET_warming - lambdaETref
DeltaET_warming_nostoma = lambdaET_warming_nostoma - lambdaETref
DeltaET_warming_VPDonly = lambdaET_warming_VPDonly - lambdaETref
DeltaET_full = lambdaET_full - lambdaETref

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(CO2, DeltaET_Rn, '-b', label=r'$\delta ET_{Rn}$')
ax1.plot(CO2, DeltaET_warming, '-r', label=r'$\delta ET_{warming}$')
ax1.plot(CO2, DeltaET_warming_nostoma, '-g', label=r'$\delta ET_{warming no stoma}$')
ax1.plot(CO2, DeltaET_warming_VPDonly, '-c', label=r'$\delta ET_{VPD only}$')
ax1.plot(CO2, DeltaET_full, '-r', linewidth=2, label=r'$\delta ET_{full}$')
plt.legend(loc='best', fontsize=8)
#set(gca,'FontSize',16)
ax1.set_xlabel('[CO_2] (ppm)')
ax1.tick_params(right=True)#,axis=ax1.yaxis)
ax1.set_ylabel(r'$\delta$ ET (W/m**2)')
ax2 = fig.add_subplot(212)
ax2.plot(CO2, 100*(Delta_globalwarming/Delta-1), '-r', label=r'$\Delta$ (#)')
ax2.plot(CO2, 100*(VPD_globalwarming/Delta_globalwarming*(Delta/VPD)-1),\
         '-b', label=r'VPD / $\Delta$ (#)')
#set(gca,'FontSize',16)
ax2.set_xlabel('[CO_2] (ppm)')
ax2.tick_params(right=True)#,axis=ax2.yaxis)
plt.legend(loc='best', fontsize=8)
# ylim([0 6])
plt.savefig('%s/temp/penman.png' % os.environ['PLOTS'])
plt.show(block=False)
