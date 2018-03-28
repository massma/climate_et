import pandas as pd
from sympy import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from shared_functions import *

pft_order = ['CSH', 'MF', 'SAV', 'ENF',\
             'EBF', 'WSA', 'DBF',  'CRO', 'GRA']

plt.close('all')
def plot_second_derivative(row, ax, zero_func):
  """
  given a fucntion for second derivative, makes a plot
  vpd on x axis wue power on y axis
  """

  # _vpd = np.linspace(dfs['5'].loc[row.pft, 'vpd'],\
  #                    dfs['95'].loc[row.pft, 'vpd'])
  _vpd = np.linspace(10., 4000.)
  _w_power = [zero_func(vpd, row.g1) for vpd in _vpd]
  print(np.nanmin(_w_power))
  ax.plot(_vpd, _w_power,\
          label="%s: g1: %5.1f"\
          % (row.pft, np.round(row.g1, 1)))
  return

init_printing()
# init rad
r, sigma = symbols('r sigma')
# init aero
g_a, p, t, cp, r_air, x, g_s = symbols('g_a P T c_p R_air x g_s')
# init plant
gamma, c_s, r_star, uwue, g_1, onesix = symbols('gamma c_s R* uWUE g_1 1.6')

vpd, delta = symbols('VPD Delta')

# powers
wue_power = symbols('n')

et_true = (delta*r\
           +g_a*p/t\
           *(cp*vpd/r_air\
             -gamma*c_s*vpd**wue_power\
             /(r_star*onesix*sigma*uwue*(1+g_1/vpd**(1/2)))))\
             /(delta + gamma)
et_pm = (delta*r + g_a*p/(t*r_air)*cp*vpd)/(delta+gamma*(1+g_a/g_s))
et = symbols('ET')
rhs = et_pm.subs(g_s, r_star*t/p*onesix*(1+ g_1/vpd**(1/2))\
                 *sigma*uwue*et/(c_s*vpd**wue_power))

# step 1
lhs = et*(delta+gamma*(1+g_a/(r_star*t/p*onesix*(1+ g_1/vpd**(1/2))\
                 *sigma*uwue*et/(c_s*vpd**wue_power))))
rhs = rhs*(delta+gamma*(1+g_a/(r_star*t/p*onesix*(1+ g_1/vpd**(1/2))\
                 *sigma*uwue*et/(c_s*vpd**wue_power))))

# step 2
lhs = simplify(lhs - g_a*gamma/(r_star*t/p*onesix*(1+ g_1/vpd**(1/2))\
                 *sigma*uwue/(c_s*vpd**wue_power)))
rhs = rhs - g_a*gamma/(r_star*t/p*onesix*(1+ g_1/vpd**(1/2))\
                 *sigma*uwue/(c_s*vpd**wue_power))

#step 3
lhs = lhs/(delta + gamma)
rhs = rhs/(delta + gamma)

#check

if simplify(rhs-et_true) == 0:
  et = et_true
else:
  raise ValueError("Algebra error somewhere")

def sgn(expr):
  """crudely returns the sign term given the d_et derivative"""
  sign = expr*t*(delta+gamma)/(p*g_a)
  neg = sign- cp/r_air
  coll = collect(neg, (c_s*gamma/(1.6*r_star*sigma*uwue*onesix)))
  return coll

d_vpd = diff(et, vpd)
# d_vpd = diff(et, e_s)+diff(et, rh)
sign = sgn(d_vpd)
d2_vpd = diff(d_vpd, vpd)
# d2_vpd = diff(d_vpd, e_s) + diff(d_vpd, rh)
et
vpd
simple_second = simplify(d2_vpd)
sign_terms = simple_second\
             *(r_star*t*vpd**6*uwue*(delta + gamma)*(vpd**(1/2)+g_1)**3)
sign_terms = simplify(sign_terms\
                      /(p*vpd**(wue_power+1/2+4)*c_s*g_a*gamma))
# sign_func = lambdify([wue_power, vpd, g_1], sign_terms)

# plot_second_derivative(sign_func, 'medlyn')
mean_df['pft'] = mean_df.index
zeros = solve(d2_vpd, wue_power)
first = zeros[0]# this solution unphysical, n = complex infinity/log(VPD)
second = zeros[1] # this solution unphysical,-0.5 <= n <= 0
third = zeros[2] # this solutions is physically reasonable, with 0.5 <= n <= 1
uwue_n_zeros = lambdify([vpd,  g_1],\
                        third)

d2_uwue = simple_second.subs(wue_power, 1/2)
zeros_uwue = solve(d2_uwue, g_1)

fig = plt.figure()
ax = fig.add_subplot(111)
for pft in pft_order:
  plot_second_derivative(mean_df.loc[pft, :], ax, uwue_n_zeros)
ax.set_xlabel('VPD (Pa)', fontsize=fontsize)
ax.set_ylabel(r'WUE VPD power ($\frac{GPP \cdot VPD^n}{ET}$)',\
              fontsize=fontsize)
ax.set_ylim([0.5, 1.0])
ax.set_xlim([0.0, 4000.])
plt.legend(loc="best")
ax.text(2500., 0.9, 'Concave Down', horizontalalignment='center',\
        verticalalignment='center', fontdict={'fontsize' : 25})

ax.text(2500., 0.55, 'Concave Up', horizontalalignment='center',\
        verticalalignment='center', fontdict={'fontsize' : 25})

plt.savefig("../../doc/paper/concave.pdf")

 
def print_range(name, func):
  """prints out limits of function"""
  print("\n****%s solution****" % name, func)
  print("limit of uwue n as g1 to infinity: ", simplify(limit(func, g_1, oo)))
  print("limit of uwue n as g1 to zero: ", simplify(limit(func, g_1, 0)))
  print("limit of uwue n as vpd to infinity: ", simplify(limit(func, vpd, oo)))
  print("limit of uwue n as vpd to zero: ", simplify(limit(func, vpd, 0)))
  return

# print_range("first", first)
# print_range("second", second)
# print_range("third", third)

et_true = (delta*r\
           +g_a*p/t\
           *(cp*vpd/r_air\
             -gamma*c_s*vpd**wue_power\
             /(r_star*onesix*sigma*uwue*(1+g_1/vpd**wue_power))))\
             /(delta + gamma)
d_vpd = diff(et_true, vpd)
d2_vpd = simplify(diff(d_vpd, vpd))
soln = solve(d2_vpd.subs((vpd**wue_power + g_1), x), x)
normed = [sol/g_1 - 1 for sol in soln]
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
n = np.linspace(0.5, 0.99)
for i, norm in enumerate(normed):
  func = lambdify([wue_power], norm)
  # ax.plot([func(_) for _ in n], n, label="sol'n %i-2" % i)
  ax.plot([(func(_)*mean_df.g1.mean())**2 for _ in n], n, label="sol'n %i" % i)
ax.set_ylabel("VPD Exponent (n)")
ax.set_xlabel(r"VPD (Pa), assuming g1=110 Pa$^{0.5}$")
ax.set_xlim([0., 4000.0])
ax.set_ylim([0.5, 1.0])
def un_normalized(normalized):
  """transform vpd**n/g1 to vpd, assuming n=0.5 and g1=100"""
  return np.round((normalized*mean_df.g1.mean())**2, 0)
def normalized(unnormalized):
  """inverse of un-normalized"""
  return np.round(np.sqrt(unnormalized)/mean_df.g1.mean(), 3)
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels(normalized(ticks))
ax2.set_xlabel(r"$\frac{VPD^n}{g_1}$")
plt.savefig('./temp.png')




et_true = (delta*r\
           +g_a*p/t\
           *(cp*vpd/r_air\
             -gamma*c_s*vpd**(1/2)\
             /(r_star*onesix*sigma*uwue*(1+g_1/vpd**wue_power))))\
             /(delta + gamma)
d_vpd = diff(et_true, vpd)
d2_vpd = simplify(diff(d_vpd, vpd))
soln = solve(d2_vpd.subs((vpd**wue_power + g_1), x), x)
normed = [sol/g_1 - 1 for sol in soln]
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
n = np.linspace(0.5, 0.99)
for i, norm in enumerate(normed):
  func = lambdify([wue_power], norm)
  # ax.plot([func(_) for _ in n], n, label="sol'n %i-2" % i)
  ax.plot([(func(_)*mean_df.g1.mean())**2 for _ in n], n, label="sol'n %i" % i)
ax.set_ylabel("VPD Exponent (n)")
ax.set_xlabel(r"VPD (Pa), assuming g1=110 Pa$^{0.5}$")
ax.set_xlim([0., 4000.0])
ax.set_ylim([0.5, 1.0])
def un_normalized(normalized):
  """transform vpd**n/g1 to vpd, assuming n=0.5 and g1=100"""
  return np.round((normalized*mean_df.g1.mean())**2, 0)
def normalized(unnormalized):
  """inverse of un-normalized"""
  return np.round(np.sqrt(unnormalized)/mean_df.g1.mean(), 3)
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels(normalized(ticks))
ax2.set_xlabel(r"$\frac{VPD^n}{g_1}$")
plt.savefig('./temp_2.png')



et_true = (delta*r\
           +g_a*p/t\
           *(cp*vpd/r_air\
             -gamma*c_s*vpd**wue_power\
             /(r_star*onesix*sigma*uwue*(1+g_1/vpd**(1/2)))))\
             /(delta + gamma)
d_vpd = diff(et_true, vpd)
d2_vpd = simplify(diff(d_vpd, vpd))
soln = solve(d2_vpd.subs((vpd**0.5 + g_1), x), x)
normed = [sol/g_1 - 1 for sol in soln]
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
n = np.linspace(0.5, 0.99)
for i, norm in enumerate(normed):
  func = lambdify([wue_power], norm)
  # ax.plot([func(_) for _ in n], n, label="sol'n %i-2" % i)
  ax.plot([(func(_)*mean_df.g1.mean())**2 for _ in n], n, label="sol'n %i" % i)
ax.set_ylabel("VPD Exponent")
ax.set_xlabel(r"VPD (Pa), assuming g1=110 Pa$^{0.5}$")
ax.set_xlim([0., 4000.0])
ax.set_ylim([0.5, 1.0])
def un_normalized(normalized):
  """transform vpd**n/g1 to vpd, assuming n=0.5 and g1=100"""
  return np.round((normalized*mean_df.g1.mean())**2, 0)
def normalized(unnormalized):
  """inverse of un-normalized"""
  return np.round(np.sqrt(unnormalized)/mean_df.g1.mean(), 3)
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels(normalized(ticks))
ax2.set_xlabel(r"$\frac{VPD^n}{g_1}$")
plt.savefig('./temp_3.png')

