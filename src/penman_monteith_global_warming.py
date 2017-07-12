import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import metcalcs as met
# penamn monteith global warming impact calculation
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

mpl.rc('font', **font)

plt.close('all')
#----------------grass set--------------------
#d = 2/3                          # zero plane displacement height
#zm = 2                           # height of wind measurements in m
#zom = 0.123                      # roughness length govern momentum transfer
#zh = 2                           # height of humidity measures
#zoh = 0.1*zom                    # roughness length hv if not given
#uz = 1                           # wind speed at height
#h =0.12                          # crop height
#k = 0.41
#rl= 70                           # single grass leaf stomatal resistance

# note that the FAO grass reference has h = 0.12, rl=70,and albedo =0.23.
#-----------------------------------------

#----------------- test set---------------------------
# unit for G_0 should be also m/s, while unit for G_1 varies with exponent of VPD (kPa)

# the mean maximum G_0 (0.35-0.4) for EBF (evergreen broadleaf forests) is 7.91*10**-3 m/s
# the mean maximum G_0 (0.35-0.4) for ENF (evergreen needleleaf forests) is 3.50*10**-3 m/s
# the mean maximum G_0 (0.35-0.4) for DBF (deciduous broadleaf forests) is 3.54*10**-3 m/s

# the mean G_1 for EBF (evergreen broadleaf forests) is about 0.88*10**-3
# the mean G_1 for ENF (evergreen needleleaf forests) is about 1.12*10**-3
# the mean G_1 for DBF (deciduous broadleaf forests) is about 1.12*10**-3

CO2  = np.arange(350.,1200.) #[ppm]
DeltaCO2 = CO2[-1]-CO2[0]
dT_dCO2  = 3.7/DeltaCO2
dRn_dCO2 = 8.5/DeltaCO2

Rn = 300. # net radiation (W/m2)
G = 0.05*Rn # soil heat flux (W/m2)
G0 = 7.91e-3 
G1 = 0.88e-3 
rhoA = 1.205 # mean air density
cP = 1004. # specific heat air
Ta = 25. # C 273.15+
RH = 70. 
# VPD = es*(1-RH/100)
k = 0.41 # vonkarmans constant
uz = 2. # wind speed at height
h = 10. # crop height
d = 2./3.*h
z0m = 0.1*h
z0h = z0m
rlmin = 20.
zmeas = 2.+h
LAI = 1. # max feasible LAI
vp_factor = 100.

ra = (np.log((zmeas-d)/z0m)*np.log((zmeas-d)/z0h))/k**2/uz
Delta = (met.vapor_pres(Ta+0.1)-met.vapor_pres(Ta))/0.1*vp_factor
es = met.vapor_pres(Ta)*vp_factor
gamma = 66. # psychrometric constant

VPD     = es*(1-RH/100)
rl      = 1./(G0+G1/np.sqrt(VPD))             # bulk stomatal resist of illuminated leaf
rs      = rl/LAI
        
lambdaETref= (Delta*(Rn-G)+rhoA*cP*VPD/ra)/(Delta+(1+gamma)*(rs/ra))
lambdaET_Rn= (Delta*(Rn+dRn_dCO2*(CO2-CO2[1])-G)+rhoA*cP*VPD/ra)/(Delta+(1+gamma)*(rs/ra))

T_globalwarming = Ta + dT_dCO2*(CO2-CO2[1])
Delta_globalwarming = (met.vapor_pres(T_globalwarming+0.1)-met.vapor_pres(T_globalwarming))/0.1*vp_factor
VPD_globalwarming   = met.vapor_pres(T_globalwarming)*(1-RH/100)*vp_factor
rs_globalwarming    = 1/(G0+G1/np.sqrt(VPD_globalwarming))/LAI

lambdaET_warming= (Delta_globalwarming*(Rn-G)+rhoA*cP*VPD_globalwarming/ra)/(Delta_globalwarming+(1+gamma)*(rs_globalwarming/ra))
lambdaET_warming_nostoma= (Delta_globalwarming*(Rn-G)+rhoA*cP*VPD_globalwarming/ra)/(Delta_globalwarming+(1+gamma)*(rs/ra))
lambdaET_warming_VPDonly= (Delta*(Rn-G)+rhoA*cP*VPD_globalwarming/ra)/(Delta+(1+gamma)*(rs_globalwarming/ra))
lambdaET_full   = (Delta_globalwarming*(Rn+dRn_dCO2*(CO2-CO2[1])-G)+rhoA*cP*VPD_globalwarming/ra)/(Delta_globalwarming+(1+gamma)*(rs_globalwarming/ra))

DeltaET_Rn = lambdaET_Rn - lambdaETref
DeltaET_warming = lambdaET_warming - lambdaETref
DeltaET_warming_nostoma = lambdaET_warming_nostoma - lambdaETref
DeltaET_warming_VPDonly = lambdaET_warming_VPDonly - lambdaETref
DeltaET_full = lambdaET_full - lambdaETref

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(CO2,DeltaET_Rn,'-b',label='$\delta ET_{Rn}$')
ax1.plot(CO2,DeltaET_warming,'-r',label='$\delta ET_{warming}$')
ax1.plot(CO2,DeltaET_warming_nostoma,'-g',label='$\delta ET_{warming no stoma}$')
ax1.plot(CO2,DeltaET_warming_VPDonly,'-c',label='$\delta ET_{VPD only}$')
ax1.plot(CO2,DeltaET_full,'-r',linewidth=2,label='$\delta ET_{full}$')
plt.legend(loc='best',fontsize=8)
#set(gca,'FontSize',16)
ax1.set_xlabel('[CO_2] (ppm)')
ax1.tick_params(right=True)#,axis=ax1.yaxis)
ax1.set_ylabel('$\delta$ ET (W/m**2)')
ax2 = fig.add_subplot(212)
ax2.plot(CO2,100*(Delta_globalwarming/Delta-1),'-r',label='$\Delta$ (#)')
ax2.plot(CO2,100*(VPD_globalwarming/Delta_globalwarming *(Delta/VPD)-1),'-b',label='VPD / $\Delta$ (#)')
#set(gca,'FontSize',16)
ax2.set_xlabel('[CO_2] (ppm)')
ax2.tick_params(right=True)#,axis=ax2.yaxis)
plt.legend(loc='best',fontsize=8)
# ylim([0 6])
plt.savefig('%s/temp/penman.png' % os.environ['PLOTS'])
plt.show(block=False)
