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
  ax.plot(t_a, np.ones(t_a.shape)*_df.scaling.mean(), 'k--',\
          linewidth=0.5, label='Term 1 mean')
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    g_a = _df.g_a.quantile(q=percentile/100.)
    scale = d_calc.scaling(mean__df, t_a, g_a)
    ax.plot(t_a, scale, label='$g_a$ = %5.3f (%dth percentile)'\
            % (g_a, int(percentile)))
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Term 1 ($\frac{g_a \; P}{T(\Delta + \gamma)}$)')
  ax.set_title('PFT = %s' % _df.pft.iloc[0])
  plt.legend(loc='best', fontsize=8)
  if savefig:
    plt.savefig('../../doc/paper/fig04.pdf')
  return

def scaled_mean(_df):
  """just scales g_a by mean"""
  _df['scaled_g_a'] = _df.g_a/df.g_a.mean()
  return _df

def plot_scaled_scaling(_df, ax):
  """makes idealized plots of scaling as a function of g_a and T"""
  t_a = np.linspace(_df.t_a.quantile(q=0.05), _df.t_a.quantile(q=0.95))
  _df = _df.groupby('pft').apply(scaled_mean)
  for percentile in [5., 25., 50., 75., 95.]:
    g_a = _df.scaled_g_a.quantile(q=percentile/100.)
    scale = d_calc.scaling(mean_df.mean(), t_a, g_a)
    ax.plot(t_a, scale,\
            label=r'$\frac{g_a}{\overline{g_a}}$ = %5.3f (%dth percentile)'\
            % (g_a, int(percentile)))
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Normalized Term 1 '\
                r'($\frac{g_a \; P}{\overline{g_a} \; T(\Delta + \gamma)}$)')
  plt.legend(loc='best', fontsize=8)
  return

def plot_mean_pft(_df, ax, savefig=False):
  """makes idealized plots of scaling as a function of pft"""

  # pfts = ['GRA', 'ENF', 'CRO', 'DBF', 'CSH']
  pfts = ['DBF', 'ENF', 'CSH', 'CRO', 'GRA']
  for pft in pfts:

    idx = (_df.pft == pft)
    t_a = np.linspace(_df.t_a.loc[idx].quantile(q=0.05),\
                      _df.t_a.loc[idx].quantile(q=0.95))
    g_a = _df.g_a.loc[idx].mean()
    scale = d_calc.scaling(mean_df.loc[pft, :], t_a, g_a)
    ax.plot(t_a, scale, label='PFT = %s, $\overline{g_a}$ = %5.3f'\
            % (pft, g_a))
    ptiles = np.array([_df.t_a.loc[idx].quantile(q=_p/100.)\
                       for _p in [25., 50., 75.]])
    ax.plot(ptiles, d_calc.scaling(mean_df.loc[pft, :], ptiles, g_a), 'k*')
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Term 1 ($\frac{g_a \; P}{T(\Delta + \gamma)}$)')
  plt.legend(loc='best', fontsize=8)
  return

def scaling_wrapper(_df):
  """wrapper that groups by pft and does scaling plot"""
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*2)
  ax = fig.add_subplot(2, 1, 1)
  plot_scaled_scaling(_df, ax)
  ax = fig.add_subplot(2, 1, 2)
  plot_mean_pft(_df, ax)
  plt.tight_layout()
  plt.savefig('../../doc/paper/fig04.pdf')
  return


# plot_scaling(df)
# os.system('convert -append %s/climate_et/scaling/*.pdf '\
#           '../../doc/paper/fig04.pdf'\
#           % (os.environ['PLOTS']))
scaling_wrapper(df)
