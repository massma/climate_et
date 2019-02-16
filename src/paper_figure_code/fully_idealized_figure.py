#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *
import matplotlib.lines as mlines

# mean_df defined her at module-level
linewidth=2.5
dashedlinewidth=3.0
markersize=8

et_string = r'$\frac{\partial \; ET}{\partial \, VPD}$ '\
                  'W m$^{-2}$ Pa$^{-1}$'

vpd = np.linspace(df.vpd.quantile(q=0.05), df.vpd.quantile(q=0.95))
median_row = df.quantile(q=0.5)
linestyles = ['-','--','-.']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:3]
ptiles = [15.0, 50.0, 85.0]
pdesc = ['min', 'median', 'max']

paren_string = r'$\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_a }{1.6 \; R\; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{VPD}}'\
               r'{2 (g_1 + \sqrt{VPD})^2}\right)$'

g1_medlyn = [2.0, 4.0, 6.0]
g1s = list(map(lambda x: x*np.sqrt(1000.0), g1_medlyn))
uwue_zhou = [6.99, 9.52, 12.05]
uwues = list(map(d_io.uWUE_converter, uwue_zhou))
ts = [10.0, 20.0, 30.0]
gas = [0.015, 0.035, 0.055]
p_default =


def get_scaling_quantiles(_df):
  out = d_calc.scaling(_df)
  return list(map((lambda p: out.quantile(q=p)), [0.15, 0.5, 0.85]))

def idealized_scaling(t, g):
  _df = {'t_a_k': t + 273.15, 'p_a': 100000.0, 'gamma': 64.5}
  return d_calc.scaling(_df, t_a=t, g_a=g)

def get_scaling_quantiles_idealized():
  x = []
  for t, g_a in zip(ts[::-1], gas):
    x.append(idealized_scaling(t, g_a))
  return x

def idealized_sign(g1, uwue):
  _df = {'gamma': 64.6, 'r_moist' : 288.0, "c_a"  : 400.0}
  return d_calc.sign(_df, g1=g1, uwue=uwue, vpd=vpd)


def get_sign_range():
  return np.array([idealized_sign(g1, uwue=uwue)
                   for g1 in g1s for uwue in uwues])

def add_plots(_df, ax):
  """takes df and plots both halves of product in term 2"""
  signs = get_sign_range()
  for scaling, scaling_ptile, c in zip(get_scaling_quantiles_idealized(),
                                       ptiles,
                                       colors):
    for sign_f, sign_ptile, l in zip([np.min, np.median, np.max],
                                     ptiles,
                                     linestyles):
      p = ax.plot(vpd, scaling*sign_f(signs, axis=0), color=c, linestyle=l)
  return

def add_sign_plots(_df, ax):
  """takes df and plots both halves of product in term 2"""
  signs = get_sign_range()
  scaling = get_scaling_quantiles(df)[1]
  for g1, scaling_ptile, c in zip(g1s, ptiles, colors):
    for uwue, sign_ptile, l in zip(uwues, ptiles, linestyles):
      print("g1: %5.5f" % g1)
      print("uwue: %5.5f" % uwue)
      p = ax.plot(vpd, scaling*idealized_sign(g1, uwue), color=c, linestyle=l)
  return

def add_scaling_plots(_df, ax):
  """takes df and plots both halves of product in term 2"""
  sign = np.median(get_sign_range(), axis=0)
  for g, scaling_ptile, c in zip(gas, ptiles, colors):
    for t, sign_ptile, l in zip(ts, ptiles, linestyles):
      print("g: %5.5f" % g)
      print("t: %5.5f" % t)
      p = ax.plot(vpd, idealized_scaling(t, g)*sign, color=c, linestyle=l)
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*1.35)
ax = fig.add_subplot(1, 1, 1)
add_plots(df, ax)
ax.set_xlabel('VPD (Pa)', fontsize=single_ax_fontsize)
ax.set_ylabel(et_string, fontsize=single_ax_fontsize)
xlims = ax.get_xlim()
ax.plot(xlims, [0., 0.], 'k--', linewidth=linewidth)
ax.set_xlim(xlims)
color_handles = [mlines.Line2D([], [], color=c, label=r'%s: ($\frac{2 \; g_a \; P}{T(\Delta + \gamma)}$)' % p) for c,p in zip(colors, pdesc)]
style_handles = [mlines.Line2D([], [], color='k', linestyle=l, label=(p + ": " + paren_string)) for l,p in zip(linestyles, pdesc)]

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*1.35)
ax = fig.add_subplot(1, 1, 1)
add_sign_plots(df, ax)
ax.set_xlabel('VPD (Pa)', fontsize=single_ax_fontsize)
# ax.set_ylabel(paren_string, fontsize=single_ax_fontsize)
ax.set_ylabel(et_string, fontsize=single_ax_fontsize)
xlims = ax.get_xlim()
ax.plot(xlims, [0., 0.], 'k--', linewidth=linewidth)
ax.set_xlim(xlims)
color_handles = [mlines.Line2D([], [], color=c, label=r'%s: g1' % p) for c,p in zip(colors, pdesc)]
style_handles = [mlines.Line2D([], [], color='k', linestyle=l, label=(p + ": uWUE")) for l,p in zip(linestyles, pdesc)]
ax.legend(handles = (color_handles + style_handles), loc='best', fontsize=single_ax_fontsize-4)
plt.tight_layout()
plt.savefig('../../doc/paper/fully_idealized_sign.pdf')

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*1.35)
ax = fig.add_subplot(1, 1, 1)
add_scaling_plots(df, ax)
ax.set_xlabel('VPD (Pa)', fontsize=single_ax_fontsize)
# ax.set_ylabel(paren_string, fontsize=single_ax_fontsize)
ax.set_ylabel(et_string, fontsize=single_ax_fontsize)
xlims = ax.get_xlim()
ax.plot(xlims, [0., 0.], 'k--', linewidth=linewidth)
ax.set_xlim(xlims)
color_handles = [mlines.Line2D([], [], color=c, label=r'%s: g_a' % p) for c,p in zip(colors, pdesc)]
style_handles = [mlines.Line2D([], [], color='k', linestyle=l, label=(p + ": T")) for l,p in zip(linestyles, pdesc)]
ax.legend(handles = (color_handles + style_handles), loc='best', fontsize=single_ax_fontsize-4)
plt.tight_layout()
plt.savefig('../../doc/paper/fully_idealized_scaling.pdf')

for p in [0.15, 0.5, 0.85]:
  print("\n*****Percentile: %5.2f" % p)
  for var in ["p_a", "g_a","t_a","gamma","delta"]:
    print("%s: %f" % (var, df[var].quantile(q=p)))
