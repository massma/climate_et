#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *

# mean_df defined her at module-level


linewidth=2.5
dashedlinewidth=3.0
markersize=8

et_string = r'$\frac{\partial \; ET}{\partial \, VPD}$ '\
                  'W m$^{-2}$ Pa$^{-1}$'

def get_scaling_quantiles(_df):
  out = d_calc.scaling(_df)
  return list(map((lambda p: out.quantile(q=p)), [0.25, 0.5, 0.75]))

def get_sign_quantiles(_df):
  def vpd_func(p, vpd):
    out = d_calc.sign(_df, vpd=vpd)
    return out.quantile(q=p)
  return list(map(lambda p: lambda vs: list(map(lambda v: vpd_func(p, v), vs)), [0.25, 0.5, 0.75]))

vpd = np.linspace(df.vpd.quantile(q=0.05), df.vpd.quantile(q=0.95))
median_row = df.quantile(q=0.5)

# for g1 in map(lambda x: x*np.sqrt(1000.0), [2.0, 4.0, 6.0]):
#   for uwue in map(d_io.uWUE_converter, [6.99, 9.52, 12.05]):

def add_plots(_df, ax):
  """takes df and plots both halves of product in term 2"""
  for scaling, scaling_ptile in zip(get_scaling_quantiles(df), [25., 50., 75.]):
    for vpd_func, sign_ptile in zip(get_sign_quantiles(df), [25., 50., 75.]):
      print(vpd_func(vpd))
      p = ax.plot(vpd, scaling*np.array(vpd_func(vpd)),\
                  label="Scaling=%3.0f Percentile, Sign= %3.0f Percentile"\
                  % (scaling_ptile, sign_ptile))
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
ax.legend(loc='best', fontsize=single_ax_fontsize-4)
plt.tight_layout()
plt.savefig('../../doc/paper/fully_idealized.pdf')
