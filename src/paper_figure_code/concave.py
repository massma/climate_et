import pandas as pd
from sympy import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from shared_functions import *

plt.close('all')

init_printing()
# init rad
r, sigma = symbols('r sigma')
# init aero
g_a, p, t, cp, r_air, x, g_s = symbols('g_a P T c_p R_air x g_s')
# init plant
gamma, c_s, r_star, uwue, g_1, onesix = symbols('gamma c_s R* uWUE g_1 1.6')

vpd, delta = symbols('VPD Delta')

def un_normalized(normalized):
  """transform vpd**n/g1 to vpd, assuming n=0.5 and g1=100"""
  return np.round((normalized*mean_df.g1.mean())**2, 0)
def normalized(unnormalized):
  """inverse of un-normalized"""
  return np.round(np.sqrt(unnormalized)/mean_df.g1.mean(), 3)

n, m = symbols('n m')
et = (delta*r\
           +g_a*p/t\
           *(cp*vpd/r_air\
             -gamma*c_s*vpd**n\
             /(r_star*onesix*sigma*uwue*(1+g_1/vpd**m))))\
             /(delta + gamma)
d_vpd = diff(et, vpd)
d2_vpd = simplify(diff(d_vpd, vpd))
soln = solve(d2_vpd.subs((vpd**m + g_1), x), x)
soln = [sol/g_1 - 1 for sol in soln] # equals vpd**m/g1
print("\n")
for sol in soln:
  pprint(sol, use_unicode=True)

fh = open("../../doc/paper/d2_solutions.tex", "w")
fh.write("\\begin{equation}\n")
fh.write("\\frac{\\partial^2 \; ET}{\\partial \\; VPD^2} = 0 "\
         "\\quad \\forall \\quad")
fh.write("%s%s%s" % ("\\frac{VPD^m}{g1} = ", latex(soln[0]), "\n"))
fh.write("\\end{equation}\n")
fh.close()

### test
substitutes = [{n : 1/2, m : m}, {n : n, m : 1/2}, {n : m, m : m}]
solve_vars = [m, n, m]
labels = ["n=1/2, m varying", "n varying, m=1/2",\
          "n=m, both vary"]


def gen_func(substitute, solve_var, soln):
  """generates a lambdified function for plotting,
  given a starting function and a subsititute function"""
  subbed = [lambdify([solve_var], sol.subs(substitute)) for sol in soln]
  return subbed

def plot_curve(ax, subbed, label):
  """plots the curve, given normlaized func"""
  n = np.linspace(0.5, 0.9999999) # singlulatity at 1
  ax.plot([(subbed[0](_)*mean_df.g1.mean())**2 for _ in n], n, label=label)
  return ax


plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
for sub, solve_var, label in zip(substitutes, solve_vars, labels):
  subbed = gen_func(sub, solve_var, soln)
  ax = plot_curve(ax, subbed, label)
ax.set_ylabel(r"VPD Exponent (m: $\frac{g_1}{VPD^m}$, "\
              r"n: $\frac{GPP}{ET}VPD^n$)", fontsize=single_ax_fontsize)
ax.set_xlabel(r"VPD (Pa), assuming g1=110 Pa$^{0.5}$ and m=1/2",\
              fontsize=single_ax_fontsize)
ax.set_xlim([0., 4000.0])
ax.set_ylim([0.5, 1.0])
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels(normalized(ticks))
# ax2.set_xlabel(r"$\frac{VPD^m}{g_1}$         ", fontsize=single_ax_fontsize)
ax2.set_xlabel(r"VPD$^m$/g$_1$", fontsize=single_ax_fontsize)
ax.text(1050., 0.95, 'Concave Down', horizontalalignment='center',\
        verticalalignment='top', fontdict={'fontsize' : 22})
ax.text(2750., 0.53, 'Concave Up', horizontalalignment='center',\
        verticalalignment='bottom', fontdict={'fontsize' : 22})
ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.6),
          fontsize=single_ax_fontsize-4)
# plt.tight_layout()
plt.savefig("../../doc/paper/concave.pdf")



