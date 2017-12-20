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


def plot_second_derivative(function, name):
  """
  given a fucntion for second derivative, makes a plot
  vpd on x axis wue power on y axis
  """
  w_power = np.linspace(0.5, 1.0)
  plt.close('all')
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*3)
  fig.set_figwidth(fig.get_figwidth()*2)
  axs = [fig.add_subplot(3,2,i+1) for i in range(5)]
  for index, ax in zip(mean_df.index, axs):
    row = dfs['mean'].loc[index, :]
    _vpd = np.linspace(dfs['5'].loc[index, 'vpd'],\
                      dfs['95'].loc[index, 'vpd'])
    _vpd, _w_power = np.meshgrid(_vpd, w_power)
    # note should probably have g1 units move with g1 power
    zvar = function(_w_power, _vpd, row.g1)
    vmax = 1000.0#1.0*zvar.std()+zvar.mean()
    color = ax.pcolormesh(_vpd, _w_power, zvar, cmap='RdBu',\
                          vmin=-vmax, vmax=vmax)
    ax.set_xlabel('VPD')
    ax.set_ylabel('WUE VPD power')
    ax.set_title('PFT: %s' % (index))
    cbar = plt.colorbar(color, ax=ax)

  cbar.set_label('Second Derivative')
  plt.tight_layout()
  plt.savefig('%s/climate_et/analytics/%s.png'\
              % (os.environ['PLOTS'], name))
  fig.clf()
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

plot_second_derivative(sign_func, 'medlyn')

g_1b, e_s = symbols('g_1b e_s')
et_bb = (delta*r\
         +g_a*p/t\
         *(cp*vpd/r_air\
           -gamma*c_s*vpd**wue_power\
           /(r_star*uwue*g_1b*(1-vpd/e_s))))\
           /(delta + gamma)

d_vpd_bb = diff(et_bb, vpd)
sign_bb = sgn(d_vpd_bb)
d2_vpd_bb = diff(d_vpd_bb, vpd)
simple = simplify(d2_vpd_bb)
sign_terms_bb = simple*r_star*t*vpd**3*g_1b*uwue*(delta+gamma)*(vpd-e_s)**3\
                /(p*c_s*e_s*g_a*gamma)

sign_func_bb = lambdify([wue_power, vpd, e_s], sign_terms_bb)
name = 'bb'
function = sign_func_bb
w_power = np.linspace(0.5, 1.0)
plt.close('all')
fig = plt.figure()
fig.set_figheight(fig.get_figheight()*3)
fig.set_figwidth(fig.get_figwidth()*2)
axs = [fig.add_subplot(3,1,i+1) for i in range(3)]
for key, ax in zip(['5', 'mean', '95'], axs):
  e_s = dfs[key].e_s.mean()
  _vpd = np.linspace(dfs['5'].vpd.min(),\
                     np.min([dfs['95'].vpd.max(), e_s]))
  _vpd, _w_power = np.meshgrid(_vpd, w_power)

  zvar = function(_w_power, _vpd, e_s)
  vmax = 1.e6#1.0*zvar.std()+zvar.mean()
  color = ax.pcolormesh(_vpd, _w_power, zvar, cmap='RdBu',\
                        vmin=-vmax, vmax=vmax)
  ax.set_xlabel('VPD')
  ax.set_ylabel('WUE VPD power')
  ax.set_title('e_s percentile: %s, value: %f' % (key, e_s))
  cbar = plt.colorbar(color, ax=ax)

cbar.set_label('Second Derivative')
plt.tight_layout()
plt.savefig('%s/climate_et/analytics/%s.png'\
            % (os.environ['PLOTS'], name))
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
