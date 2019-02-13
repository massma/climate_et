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
ptiles = [25.0, 50.0, 75.0]
pdesc = ['min', 'median', 'max']

paren_string = r'$\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_a }{1.6 \; R\; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{VPD}}'\
               r'{2 (g_1 + \sqrt{VPD})^2}\right)$'

g1s = list(map(lambda x: x*np.sqrt(1000.0), [2.0, 4.0, 6.0]))
uwues = list(map(d_io.uWUE_converter, [6.99, 9.52, 12.05]))
def get_scaling_quantiles(_df):
  out = d_calc.scaling(_df)
  return list(map((lambda p: out.quantile(q=p)), [0.25, 0.5, 0.75]))


def get_sign_range():
    return np.array([d_calc.sign(median_row, g1=g1, uwue=uwue, vpd=vpd)
                     for g1 in g1s for uwue in uwues])

def add_plots(_df, ax):
  """takes df and plots both halves of product in term 2"""
  signs = get_sign_range()
  for scaling, scaling_ptile, c in zip(get_scaling_quantiles(df),
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
      p = ax.plot(vpd, scaling*d_calc.sign(median_row, g1=g1, uwue=uwue, vpd=vpd), color=c, linestyle=l)
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
