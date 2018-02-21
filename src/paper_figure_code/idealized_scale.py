#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *

### note mean_df defined at this level)

### FIGURE 4 ###
# scaling term as a function of t and g_a
for key in ['gamma', 'delta', 'rho_a', 'g_a']:
  print('%s var: %f, mean: %f, cf: %f'\
        % (key, df[key].std(), df[key].mean(), df[key].std()/df[key].mean()))
df['rho_like'] = df['p_a']/(273.15 + df['t_a'])

def plot_scaling(_df, ax, savefig=False):
  """makes idealized plots of scaling as a function of g_a and T"""
  t_a = np.linspace(_df.t_a.quantile(q=0.05), _df.t_a.quantile(q=0.95))
  # ax.plot(t_a, np.ones(t_a.shape)*_df.scaling.mean(), 'k--',\
  #         linewidth=0.5, label='Term 1 mean')
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    g_a = _df.g_a.quantile(q=percentile/100.)
    scale = d_calc.scaling(mean_df.loc[_df.pft.iloc[0]], t_a=t_a, g_a=g_a)
    ax.plot(t_a, scale, label='$g_a$ = %5.3f (%dth percentile)'\
            % (g_a, int(percentile)))
  custom_xlabel(_df, ax, 'T (C)')
  custom_ylabel(_df, ax,\
                r'($\frac{g_a \; P}{T(\Delta + \gamma)}$)')
  ax.set_title(name_dict[_df.pft.iloc[0]], fontsize=fontsize+3)
  plt.legend(loc='best', fontsize=10)
  if savefig:
    plt.savefig('../../doc/paper/idealized_scale.pdf')
  return

def scaling_wrapper(_df):
  """wrapper that groups by pft and does scaling plot"""
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*3)
  fig.set_figwidth(fig.get_figwidth()*2)
  for i, pft in enumerate(pft_order):
    print(pft)
    ax = fig.add_subplot(3, 2, i+1)
    plot_scaling(_df.loc[(_df.pft==pft), :], ax)
  plt.tight_layout()
  plt.savefig('../../doc/paper/idealized_scale.pdf')
  return

scaling_wrapper(df)
