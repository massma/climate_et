#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *
mpl.rcParams.update(small_ax_params)
### note mean_df defined at this level)

### FIGURE 4 ###
# scaling term as a function of t and g_a
for key in ['gamma', 'delta', 'rho_a', 'g_a']:
  print('%s var: %f, mean: %f, cf: %f'\
        % (key, df[key].std(), df[key].mean(), df[key].std()/df[key].mean()))
df['rho_like'] = df['p_a']/(273.15 + df['t_a'])

# set lims
ylim = [0.0, 0.6]
xlim = t_lims

def plot_scaling(_df, ax, savefig=False):
  """makes idealized plots of scaling as a function of g_a and T"""
  t_a = np.linspace(_df.t_a.quantile(q=0.05), _df.t_a.quantile(q=0.95))
  # ax.plot(t_a, np.ones(t_a.shape)*_df.scaling.mean(), 'k--',\
  #         linewidth=0.5, label='Term 1 mean')
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    g_a = _df.g_a.quantile(q=percentile/100.)
    scale = d_calc.scaling(mean_df.loc[_df.pft.iloc[0]], t_a=t_a, g_a=g_a)
    ax.plot(t_a, scale, label='$g_a$ = %5.3f (%dth p-tile)'\
            % (g_a, int(percentile)))
  custom_xlabel(_df, ax, 'T (C)', fontsize=small_ax_fontsize)
  custom_ylabel(_df, ax,\
                r'($\frac{2 \; g_a \; P}{T(\Delta + \gamma)}$)',\
                fontsize=small_ax_fontsize+4)
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  ax.set_title(name_dict[_df.pft.iloc[0]], fontsize=fontsize+3)
  plt.legend(loc='best', fontsize=fontsize-2)
  if savefig:
    plt.savefig('../../doc/paper/idealized_scale.pdf')
  return


panel_wrapper(df, plot_scaling, "idealized_scale.pdf")
