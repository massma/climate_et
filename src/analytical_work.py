import pandas as pd
from sympy import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np

with open('%s/changjie/diagnosed_data.pkl' % os.environ['DATA'],\
          mode='rb') as file:
  dfs = pickle.load(file)
df = dfs['full']
mean_df = dfs['mean']


# #try some anayltics with sympy
# g_1 = Symbol('g_1')
# x = Symbol('\frac{\sqrt{D}}{g_1}')
# func = g_1*(2 + x)/(2*g_1**2*(1 + x)**2)
# print(latex(series(func, x, x0=0., dir='+', n=4)))
init_printing()
# init rad
r, delta = symbols('r Delta')
# init aero
g_a, p, t, cp, vpd, r_air = symbols('g_a P T c_p VPD R_air')
# init plant
gamma, c_s, r_star, uwue, g_1 = symbols('gamma c_s R* uWUE g_1')

# powers
wue_power, g1_power = symbols('k n')

et = (delta*r\
      +g_a*p/t\
       *(cp*vpd/r_air\
         -gamma*c_s*vpd**wue_power\
          /(r_star*1.6*uwue*(1+g_1/vpd**g1_power))))\
      /(delta + gamma)

d_vpd = diff(et, vpd)
d2_vpd = diff(d_vpd, vpd)
et
vpd
simple_second = simplify(d2_vpd)
sign_terms = simple_second\
             *(r_star*t*vpd**6*uwue*(delta + gamma)*(vpd**g1_power+g_1)**3)
sign_terms = simplify(sign_terms\
                      /(p*vpd**(wue_power+g1_power+4)*c_s*g_a*gamma))
sign_func = lambdify([wue_power, g1_power, g_1, vpd], sign_terms)

w_power = np.linspace(0.5, 1.0)
g_power = np.linspace(0.5, 1.0)
w_power, g_power = np.meshgrid(w_power, g_power)
plt.close('all')
fig = plt.figure()
fig.set_figheight(fig.get_figheight()*3)
for index in mean_df.index:
  axs = [fig.add_subplot(3,1,i+1) for i in range(3)]
  for ax, key in zip(axs, ['5', 'mean', '95']):
    row = dfs[key].loc[index, :]
    # note should probably have g1 units move with g1 power
    g1_corrected = (row.g1**2)**g_power
    zvar = sign_func(w_power, g_power, g1_corrected, row.vpd)
    vmax = 1000.0#1.0*zvar.std()+zvar.mean()
    color = ax.pcolormesh(w_power, g_power, zvar, cmap='RdBu',\
                          vmin=-vmax, vmax=vmax)
    ax.set_xlabel('WUE VPD power')
    ax.set_ylabel('Stomatal VPD power')
    ax.set_title('PFT: %s, %s percentile VPD' % (index, key))
    plt.colorbar(color, ax=ax)
  plt.savefig('%s/climate_et/analytics/%s.png' % (os.environ['PLOTS'], index))
  fig.clf()

# #belos is doing split of d_s and vpd
# # vpd = (1. - rh)*e_s
# t_a = Symbol('T')
# rh = Symbol('RH')

# e_s = Symbol('e_s')# 610.8*exp((17.27*t_a)/(237.3 + t_a))
# vpd = Symbol('D_s')

# e_s = Function('e_s')(t_a)
# rh_func = (1 - vpd/e_s)
# e_s_func = vpd/(1 - rh)
# delta_es = Derivative(e_s, t_a)

# print(latex(simplify(delta_es)))
# et = (delta_es*r + g_a*rho*(c_p*(1-rh)*e_s/r_air\
#                          -gamma*c_s*sqrt((1 - rh)*e_s)\
#                          /(r_star*1.6*sigma*uwue*(1+g_1/sqrt((1-rh)*e_s)))))\
#                          /(delta_es + gamma)

# print(diff(rh_func, vpd))
# print(diff(et, rh))
# test = diff(et, e_s)

# deriv = simplify(diff(et, rh)*diff(rh_func, vpd)\
#                  +diff(et, e_s)*diff(e_s_func, vpd))
# print('full ', deriv)
# delta = Symbol('Delta')
# subbed = deriv.subs(delta_es, delta)
# subbed_2 = subbed.subs(e_s*(rh-1), -vpd)
# subbed_3 = simplify(subbed_2)
# subbed_3

# print('latex:\n', latex(subbed_3))

# e_s = Function('e_s')(t_a)
# delta_es = Derivative(e_s, t_a)
# term1 = g_a*rho/(gamma + delta_es)
# test2 = Derivative(term1, e_s)

# test2


# #below was originally in shared_functions, computes series:
# #try some analytis with sympy
# g1 = Symbol('g_1')
# x = Symbol('\frac{\sqrt{D}}{g_1}')
# func = g1*(2 + x)/(2*g1**2*(1 + x)**2)
# print(latex(series(func, x, x0=0., dir='+', n=4)))
# # x = Symbol('\sqrt{D}')
# # func  = (2*g1 + x)/(2*(g1 + x)**2)
# # print(latex(series(func, x, x0=0., dir='+')))
