#! ~/edward/bin/python
"""
This is a module to test if vpd solver is working.
"""
import codebase.penman_monteith as pm
from scipy.optimize import fsolve
import numpy as np
import metcalcs as met
import importlib

importlib.reload(pm)

atmos = {}
atmos['r_n'] = 300. #W/m2
atmos['rho_a'] = 1.205 #density kg/m3
atmos['t_a'] = 25. # C
atmos['rh'] = 70. # rel humdidity

atmos['u_z'] = 2. #wind speed at meas. hiehgt (m/s)
atmos['vpd'] = met.vapor_pres(atmos['t_a'])*(1.-atmos['rh']/100.)*100.
atmos['co2'] = 400. # ppm

canopy = {}
canopy['pft'] = 'DBF'
canopy['height'] = 10. # plant heigh m
canopy['lai'] = 1. # leaf area index pierre says 1 max feasible

et = pm.medlyn_penman_monteith(atmos, canopy)

print('Medlyn r_l: %f, Oren r_l: %f' % \
      (atmos['r_s'], pm.oren_r_l(atmos['vpd'], canopy['pft'])))
atmos.pop('r_s')
et_oren = pm.penman_monteith(atmos, canopy)
print('Medlyn et: %f, Oren et: %f' % \
      (et, et_oren))

# _et = pm.optimizer_wrapper(100.,atmos,canopy)
