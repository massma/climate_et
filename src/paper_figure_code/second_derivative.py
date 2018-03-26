import pandas as pd
from sympy import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from shared_functions import *

def plot_second_derivative(row, ax, d2_func):
  """
  given a fucntion for second derivative, makes a plot
  vpd on x axis wue power on y axis
  """
  if row.shape[0] == 1:
    row = row.iloc[0]
  else:
    raise ValueError("Something other than row passed to second deriv.",\
                     row.shape)
  w_power = np.linspace(0.5, 1.0)
  _vpd = np.linspace(dfs['5'].loc[row.pft, 'vpd'],\
                     dfs['95'].loc[row.pft, 'vpd'])
  _vpd, _w_power = np.meshgrid(_vpd, w_power)
  # note should probably have g1 units move with g1 power
  #zvar = function(_w_power, _vpd, row.g1)
  zvar = d2_func(_w_power, row.g_a, row.p_a, row.t_a_k, d_calc.CP, _vpd,\
                 row.r_moist, row.gamma, row.delta, row.r_n,\
                 row.c_a, d_calc.R_STAR, row.uwue, row.g1, 1.6)

  vmax = 1.e-5 #1000.0#1.0*zvar.std()+zvar.mean()
  # vmax = 1000.0#1.0*zvar.std()+zvar.mean()
  color = ax.contourf(_vpd, _w_power, zvar, cmap='Greys',\
                        vmin=-vmax, vmax=vmax)
  custom_xlabel(row.pft, ax, 'VPD (Pa)', fontsize=small_ax_fontsize)
  custom_ylabel(row.pft, ax, 'WUE VPD power', fontsize=small_ax_fontsize)

  # ax.set_xlabel('VPD')
  # ax.set_ylabel('WUE VPD power')
  ax.set_title(name_dict[row.pft],\
               fontsize=fontsize+3)

  return

init_printing()
# init rad
r, delta = symbols('r Delta')
# init aero
g_a, p, t, cp, vpd, r_air = symbols('g_a P T c_p VPD R_air')
# init plant
gamma, c_s, r_star, uwue, g_1, onesix = symbols('gamma c_s R* uWUE g_1 1.6')

# powers
wue_power = symbols('n')

et = (delta*r\
      +g_a*p/t\
       *(cp*vpd/r_air\
         -gamma*c_s*vpd**wue_power\
          /(r_star*onesix*uwue*(1+g_1/vpd**(1/2)))))\
      /(delta + gamma)

def sgn(expr):
  """crudely returns the sign term given the d_et derivative"""
  sign = expr*t*(delta+gamma)/(p*g_a)
  neg = sign- cp/r_air
  coll = collect(neg, (c_s*gamma/(1.6*r_star*uwue*onesix)))
  return coll

d_vpd = diff(et, vpd)
sign = sgn(d_vpd)
d2_vpd = diff(d_vpd, vpd)
et
vpd
simple_second = simplify(d2_vpd)
sign_terms = simple_second\
             *(r_star*t*vpd**6*uwue*(delta + gamma)*(vpd**(1/2)+g_1)**3)
sign_terms = simplify(sign_terms\
                      /(p*vpd**(wue_power+1/2+4)*c_s*g_a*gamma))
sign_func = lambdify([wue_power, vpd, g_1], sign_terms)

# plot_second_derivative(sign_func, 'medlyn')
mean_df['pft'] = mean_df.index
df_function = lambdify([wue_power, g_a, p, t, cp, vpd, r_air, gamma,\
                        delta, r, c_s, r_star, uwue, g_1, onesix],\
                       d2_vpd)
panel_wrapper(mean_df, plot_second_derivative, "concave.pdf",\
              args=(df_function,))

