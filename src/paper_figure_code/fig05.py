#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *

# mean_df defined her at module-level

fontsize=14
linewidth=2.5
dashedlinewidth=3.0
markersize=8
# below is whether to plot taylor series approximation
#below is to make slight changes for the talk
PLOT_TALK = False

PLOT_PET = False

name_dict = {'CRO': 'crops',\
             'DBF': 'deciduous forest',
             'ENF': 'evergreen forest',
             'GRA': 'grass',
             'CSH': 'shrub (closed)'}

if PLOT_TALK:
  paren_string = r'$\left(\frac{ c_p}{R_{air}} '\
                 r'- \frac{\gamma c_s }{1.6 \; R\; uWUE  }'\
                 r'\left( \frac{2 g_1 + \sqrt{D}}'\
                 r'{2 (g_1 + \sqrt{D})^2}\right)\right)$'
else:
  paren_string = r'(Term 2 - Term 3) $\left(\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_s }{\sigma \; 1.6 \; R\; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{D}}'\
               r'{2 (g_1 + \sqrt{D})^2}\right)\right)$'



I = 5
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
                    label="%s: $uWUE\cdot\overline{\sigma}$=%4.2f,f g1=%4.1f"\
                    % (_df.pft.iloc[0],\
                       uwue,  _df.g1.iloc[0]))

  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(vpd/1.0e3, np.ones(vpd.shape)*I, linewidth=linewidth)
  axs[1].plot(ptiles/1.0e3,\
              np.ones(ptiles.shape)*I, 'k*', markersize=markersize)
  vpd_crit = et_min_vpd(mean_row, uwue=uwue)/1.0e3
  axs[1].plot([vpd_crit], [I], 'k|', markersize=15, mew=3)
  axs[-1].plot(np.ones(vpd.shape)*I, vpd)
  axs[-1].plot(np.ones(ptiles.shape)*I, ptiles, 'k*')

  uwue = np.linspace(_df.uwue.quantile(q=0.05), _df.uwue.quantile(q=0.95))
  vpd = mean_row.vpd
  et_min = et_min_vpd(mean_row, uwue=uwue)
  axs[2].plot(uwue/mean_row.uwue_zhou, et_min,\
              label=r"PFT = %s, uWUE$\cdot\overline\sigma$=%4.2f, g1=%4.1f"\
              % (_df.pft.iloc[0],\
                 mean_row.uwue, _df.g1.iloc[0]))
  ptiles = np.array([_df.uwue.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[3].plot(uwue/mean_row.uwue_zhou, np.ones(uwue.shape)*I)
  axs[3].plot(ptiles/mean_row.uwue_zhou, np.ones(ptiles.shape)*I, 'k*')
  I -= 1
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*2.5)
_axs = [fig.add_subplot(2, 1, i+1) for i in range(2)]
axs = []
for _ax in _axs:
  divider = make_axes_locatable(_ax)
  axs.append(_ax)
  axs.append(divider.append_axes("bottom", size="20%", pad=0.0, sharex=_ax))
axs.append(divider.append_axes("left", size="20%", pad=0.0, sharey=axs[2]))
I = 5

for pft in ['CRO', 'DBF', 'GRA', 'ENF', 'CSH']:
  _df = df.loc[df.pft == pft, :]
  pft_leaf(_df, axs)

axs[1].set_xlabel('VPD (kPa)', fontsize=fontsize)
axs[3].set_xlabel('$\sigma$')
axs[0].set_ylabel(paren_string, fontsize=fontsize)
axs[0].plot(axs[0].get_xlim(), [0., 0.], 'k--', linewidth=linewidth)
if PLOT_PET:
  axs[0].plot(axs[0].get_xlim(), [d_calc.CP/287.0, d_calc.CP/287.0],\
              'm--', linewidth=dashedlinewidth, label='PET')
plt.setp(axs[2].get_yticklabels(), visible=False)
axs[-1].set_ylabel(r'VPD')#$_{ETmin}$')
axs[2].set_title(r'D$_{crit}$'\
                 r' (where $\frac{\partial \; ET}{\partial \; D} = 0$)')
axs[2].set_ylim([0., 4000.])
axs[1].set_ylim([0.5,5.5])
axs[3].set_ylim([0.5,5.5])
axs[1].get_yaxis().set_visible(False)
axs[3].get_yaxis().set_visible(False)
axs[-1].set_xlim([0.5,5.5])
axs[-1].get_xaxis().set_visible(False)

axs[2].text(0.35, 250., r'$\frac{\partial \; ET}{\partial \; D} < 0$',\
            fontsize=15)
axs[2].text(1.4, 3500., r'$\frac{\partial \; ET}{\partial \; D} > 0$',\
            fontsize=15)

if not PLOT_PET:
  axs[0].set_ylim((-1.8058955891452384, 0.87408219563044132))
for ax in axs[:1]:
  h, l = ax.get_legend_handles_labels()
  if not PLOT_PET:
    ax.legend(h, l, loc='best', fontsize=fontsize)
# plt.legend(loc='best')
plt.tight_layout()
if PLOT_TALK:
  #plt.savefig('../../doc/shared_figs/fig05_pet.pdf')
  if PLOT_PET:
    plt.savefig('../../doc/shared_figs/fig05_pet_garb.pdf')
  else:
    plt.savefig('../../doc/shared_figs/fig05_garb.pdf')
else:
  plt.savefig('../../doc/paper/fig05_garb.pdf')

