#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *

# mean_df defined her at module-level


linewidth=2.5
dashedlinewidth=3.0
markersize=8

paren_string = r'$\frac{\partial \; ET}{\partial \, VPD}$ '\
                  'W m$^{-2}$ Pa$^{-1}$'

def get_scaling_quantile(df):
  out = d_calc.scaling(df)
  return list(map((lambda p: out.quantile(q=p)), [0.25, 0.5, 0.75]))

def add_plots(_df, axs):
  """takes df and plots both halves of product in term 2"""
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  median_row = mean_df.loc[_df.pft.iloc[0], :]
  for scaling in [get_scaling_quantile(df)[1]]:
    for g1 in map(lambda x: x*np.sqrt(1000.0), [2.0, 4.0, 6.0]):
      for uwue in map(d_io.uWUE_converter, [6.99, 9.52, 12.05]):
        p = axs[0].plot(vpd, scaling*d_calc.sign(median_row, vpd=vpd, uwue=uwue, g1=g1),\
                        label="scaling=%4.1f, $uWUE$=%4.2f, g1=%4.1f"\
                        % (np.round(scaling, 3),\
                           np.round(uwue, 2),  np.round(g1, 1)))
  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(vpd, np.ones(vpd.shape), linewidth=linewidth)
  axs[1].plot(ptiles,\
              np.ones(ptiles.shape), 'k*', markersize=markersize)
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*1.35)
_axs = [fig.add_subplot(1, 1, i+1) for i in range(1)]
axs = []
for _ax in _axs:
  divider = make_axes_locatable(_ax)
  axs.append(_ax)
  axs.append(divider.append_axes("bottom", size="35%", pad=0.0, sharex=_ax))

add_plots(df, axs)

axs[1].set_xlabel('VPD (Pa)', fontsize=single_ax_fontsize)
axs[0].set_ylabel(paren_string, fontsize=single_ax_fontsize)
axs[0].plot(vpd_xlim, [0., 0.], 'k--', linewidth=linewidth)

axs[1].set_ylim([0.5,9.5])
axs[1].get_yaxis().set_visible(False)
axs[0].set_xlim(vpd_xlim)
axs[1].set_xlim(vpd_xlim)
for ax in axs[:1]:
  h, l = ax.get_legend_handles_labels()
  ax.legend(h, l, loc='best', fontsize=single_ax_fontsize-4)
plt.tight_layout()
plt.savefig('../../doc/paper/fully_idealized.pdf')
