#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *

# mean_df defined her at module-level

fontsize=10

linewidth=2.5
dashedlinewidth=3.0
markersize=8
# below is whether to plot taylor series approximation
#below is to make slight changes for the talk
PLOT_TALK = False

PLOT_PET = False

paren_string = r'$\left(\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_s }{1.6 \; R\; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{D}}'\
               r'{2 (g_1 + \sqrt{D})^2}\right)\right)$'

if PLOT_TALK:
  fontsize=14

# mask pft order from shared func to do top to bottom legnend
pft_order = ['CRO', 'DBF', 'EBF', 'ENF', 'GRA', 'CSH']

I = 6
def pft_leaf(_df, axs):
  """takes df and plots both halves of product in term 2"""
  global I
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  mean_row = mean_df.loc[_df.pft.iloc[0], :]
  uwue = mean_row.uwue
  if PLOT_TALK:
    p = axs[0].plot(vpd/1.0e3, d_calc.sign(mean_row, vpd=vpd, uwue=uwue),\
                    linewidth=linewidth,\
                    label="%s"#: uWUE=%4.2f, g1=%4.1f"\
                    % (name_dict[_df.pft.iloc[0]]))# ,\
                       # _df.uwue.iloc[0],  _df.g1.iloc[0]))
  else:
    p = axs[0].plot(vpd, d_calc.sign(mean_row, vpd=vpd, uwue=uwue),\
                    label="%s: $uWUE$=%4.2f, g1=%4.1f"\
                    % (_df.pft.iloc[0],\
                       uwue,  _df.g1.iloc[0]))

  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(vpd, np.ones(vpd.shape)*I, linewidth=linewidth)
  axs[1].plot(ptiles,\
              np.ones(ptiles.shape)*I, 'k*', markersize=markersize)
  vpd_crit = et_min_vpd(mean_row, uwue=uwue)
  axs[1].plot([vpd_crit], [I], 'k|', markersize=15, mew=3)
  # axs[-1].plot(np.ones(vpd.shape)*I, vpd)
  # axs[-1].plot(np.ones(ptiles.shape)*I, ptiles, 'k*')

  uwue = np.linspace(_df.uwue.quantile(q=0.05), _df.uwue.quantile(q=0.95))
  vpd = mean_row.vpd
  et_min = et_min_vpd(mean_row, uwue=uwue)
  # axs[2].plot(uwue/mean_row.uwue_zhou, et_min,\
  #             label=r"PFT = %s, uWUE$\cdot\overline\sigma$=%4.2f, g1=%4.1f"\
  #             % (_df.pft.iloc[0],\
  #                mean_row.uwue, _df.g1.iloc[0]))
  # ptiles = np.array([_df.uwue.quantile(q=_p/100.)\
  #                    for _p in [25., 50., 75.]])
  # axs[3].plot(uwue/mean_row.uwue_zhou, np.ones(uwue.shape)*I)
  # axs[3].plot(ptiles/mean_row.uwue_zhou, np.ones(ptiles.shape)*I, 'k*')
  I -= 1
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*1.25)
_axs = [fig.add_subplot(1, 1, i+1) for i in range(1)]
axs = []
for _ax in _axs:
  divider = make_axes_locatable(_ax)
  axs.append(_ax)
  axs.append(divider.append_axes("bottom", size="25%", pad=0.0, sharex=_ax))
# axs.append(divider.append_axes("left", size="20%", pad=0.0, sharey=axs[2]))
I = 6

for pft in pft_order:
  _df = df.loc[df.pft == pft, :]
  pft_leaf(_df, axs)

axs[1].set_xlabel('VPD', fontsize=fontsize)
#axs[3].set_xlabel('$\sigma$')
axs[0].set_ylabel(paren_string, fontsize=fontsize)
axs[0].plot(vpd_xlim, [0., 0.], 'k--', linewidth=linewidth)
if PLOT_PET:
  axs[0].plot(vpd_xlim, [d_calc.CP/287.0, d_calc.CP/287.0],\
              'm--', linewidth=dashedlinewidth, label='PET')

axs[1].set_ylim([0.5,6.5])
axs[1].get_yaxis().set_visible(False)
axs[0].set_xlim(vpd_xlim)
axs[1].set_xlim(vpd_xlim)
# if not PLOT_PET:
#   axs[0].set_ylim((-1.8058955891452384, 0.87408219563044132))
for ax in axs[:1]:
  h, l = ax.get_legend_handles_labels()
  if not PLOT_PET:
    ax.legend(h, l, loc='best', fontsize=fontsize)
# plt.legend(loc='best')
plt.tight_layout()
if PLOT_TALK:
  #plt.savefig('../../doc/shared_figs/fig04_pet.pdf')
  if PLOT_PET:
    plt.savefig('../../doc/shared_figs/fig04_pet_garb.pdf')
  else:
    plt.savefig('../../doc/shared_figs/fig04_garb.pdf')
else:
  plt.savefig('../../doc/paper/idealized_sign.pdf')

