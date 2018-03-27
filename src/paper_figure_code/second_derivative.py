import pandas as pd
from sympy import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from shared_functions import *

pft_order = ['CSH', 'MF', 'SAV', 'ENF',\
             'EBF', 'WSA', 'DBF',  'CRO', 'GRA']


def plot_second_derivative(row, ax, zero_func):
  """
  given a fucntion for second derivative, makes a plot
  vpd on x axis wue power on y axis
  """

  _vpd = np.linspace(dfs['5'].loc[row.pft, 'vpd'],\
                     dfs['95'].loc[row.pft, 'vpd'])
  def temp_f(vpd):
    return zero_func(row.g_a, row.p_a, row.t_a_k, d_calc.CP, vpd,\
                   row.r_moist, row.gamma, row.delta, row.r_n,\
                   row.c_a, d_calc.R_STAR, row.uwue, row.g1, 1.6)
  _w_power = [temp_f(vpd) for vpd in _vpd]
  ax.plot(_vpd, _w_power,\
          label="%s"\
          % row.pft)


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
zeros = solve(d2_vpd, wue_power)
uwue_n_zeros = lambdify([g_a, p, t, cp, vpd, r_air, gamma,\
                    delta, r, c_s, r_star, uwue, g_1, onesix],\
                   zeros[2])
fig = plt.figure()
ax = fig.add_subplot(111)
for pft in pft_order:
  plot_second_derivative(mean_df.loc[pft, :], ax, uwue_n_zeros)
ax.set_xlabel('VPD (Pa)', fontsize=fontsize)
ax.set_ylabel('WUE VPD power', fontsize=fontsize)
ax.set_ylim([0.5, 1.0])
plt.legend(loc="best")
ax.text(2500., 0.9, 'Concave Down', horizontalalignment='center',\
        verticalalignment='center', fontdict={'fontsize' : 25})

ax.text(2500., 0.55, 'Concave Up', horizontalalignment='center',\
        verticalalignment='center', fontdict={'fontsize' : 25})

plt.savefig("../../doc/paper/concave.pdf")
