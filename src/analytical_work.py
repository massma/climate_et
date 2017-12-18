import pandas as pd
from sympy import *

#below is just for checking of scaling
df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])

#try some anayltics with sympy
g_1 = Symbol('g_1')
x = Symbol('\frac{\sqrt{D}}{g_1}')
func = g_1*(2 + x)/(2*g_1**2*(1 + x)**2)
print(latex(series(func, x, x0=0., dir='+', n=4)))

r  = Symbol('R')

g_a = Symbol('g_a')
rho = Symbol('rho') # note that rho si not rho but P/T
c_p = Symbol('c_p')
r_air = Symbol('R_{air}')
gamma = Symbol('gamma')
c_s = Symbol('c_s')
r_star = Symbol('R*')
sigma = Symbol('sigma')
uwue = Symbol('uWUE')

t_a = Symbol('T')
rh = Symbol('RH')

e_s = Symbol('e_s')# 610.8*exp((17.27*t_a)/(237.3 + t_a))
e_s = Function('e_s')(t_a)

vpd = Symbol('D_s')
# vpd = (1. - rh)*e_s

rh_func = (1 - vpd/e_s)
e_s_func = vpd/(1 - rh)

delta_es = Derivative(e_s, t_a)

print(latex(simplify(delta_es)))
et = (delta_es*r + g_a*rho*(c_p*(1-rh)*e_s/r_air\
                         -gamma*c_s*sqrt((1 - rh)*e_s)\
                         /(r_star*1.6*sigma*uwue*(1+g_1/sqrt((1-rh)*e_s)))))\
                         /(delta_es + gamma)

print(diff(rh_func, vpd))
print(diff(et, rh))
test = diff(et, e_s)

deriv = simplify(diff(et, rh)*diff(rh_func, vpd)\
                 +diff(et, e_s)*diff(e_s_func, vpd))
print('full ', deriv)
delta = Symbol('Delta')
subbed = deriv.subs(delta_es, delta)
subbed_2 = subbed.subs(e_s*(rh-1), -vpd)
subbed_3 = simplify(subbed_2)
subbed_3

print('latex:\n', latex(subbed_3))

e_s = Function('e_s')(t_a)
delta_es = Derivative(e_s, t_a)
term1 = g_a*rho/(gamma + delta_es)
test2 = Derivative(term1, e_s)

test2


#below was originally in shared_functions, computes series:
#try some analytis with sympy
g1 = Symbol('g_1')
x = Symbol('\frac{\sqrt{D}}{g_1}')
func = g1*(2 + x)/(2*g1**2*(1 + x)**2)
print(latex(series(func, x, x0=0., dir='+', n=4)))
# x = Symbol('\sqrt{D}')
# func  = (2*g1 + x)/(2*(g1 + x)**2)
# print(latex(series(func, x, x0=0., dir='+')))
