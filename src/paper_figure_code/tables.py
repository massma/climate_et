#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""

from shared_functions import *
import scipy.optimize
import matplotlib.pyplot as plt

full_clean = True

if full_clean:
  def clean_df(_df, var='lai'):
    """remove unphysical LAI values from a df"""
    out = _df.loc[((_df[var] > 0.1) & (_df[var] < 100.)), :]
    return out

  def site_clean(_df, var='lai'):
    """this will remove some percentile of data"""
    out = _df.loc[((_df[var] < _df[var].quantile(q=0.95)) & \
                  (_df[var] > _df[var].quantile(q=0.05))), :]
    return out

  print('all removal')
  df = pd.read_pickle('%s/changjie/full_pandas_fix_scaling.pkl'\
                      % os.environ['DATA'])
  print(df.shape)
  df = site_clean(df)
  print('outlier', df.shape)
  df = clean_df(df)
  print('cleaned', df.shape)
  df = clean_df(df, var='lai_gpp')
  print('lceaned gpp', df.shape)
  print(df.shape)
  df['g_a'] = 1./df['r_a']
  df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
  df['d_et_atm'] = df['scaling']*df['vpd_atm']
  df['uwue_norm'] = df.uwue/pm.LV
  df['r_net'] = df['r_n'] - df['g_flux']
  jacobians = {'vpd' : calc.d_et,\
               'lai' : calc.d_et_d_lai,\
               'seasonal_lai' : calc.d_et_d_lai,\
               'residual_lai' : calc.d_et_d_lai,\
               'g_a' : calc.d_et_d_g_a,\
               'delta' : calc.d_et_d_delta,\
               'r_net' : calc.d_et_d_r_net}

plt.close('all')
###### table 5 ####
def frequency(_df):
  """return fraction fos amples d_et < 0"""
  return _df.d_et[_df.d_et < 0.].count() / _df.d_et.count()


pft = df.groupby('site').apply(get_pft)
mean = df.groupby('site').mean()
mean['pft'] = pft
std = df.groupby('site').std()
#jacobian = site_analysis(mean)
mean['d_et_bar'] = calc.d_et(mean)
mean['d_et_bar_std'] = np.absolute(calc.d_et(mean))*std['vpd']
mean['d_rn_bar'] = np.absolute(calc.d_et_d_r_net(mean))*std['r_net']
mean['d_et_bar_norm_rn'] = mean['d_et_bar_std']\
                           /mean['d_rn_bar']

mean_df = mean.groupby('pft').mean()
counts = df.groupby('pft').apply(frequency)

print('\n mean d_et\n', mean_df.d_et)
print('\n d_et(mean)\n', mean_df.d_et_bar)
print('\n mean d_et *std vpd\nn', mean_df['d_et_bar_std'])
print('\n mean d_et / rad\n', mean_df['d_et_bar_norm_rn'])
print('\n portion d_et < 0 \n',\
      counts)
mean_df['counts'] = counts
test = mean_df.loc[:, ['d_et', 'd_et_bar', 'd_et_bar_std',\
                       'd_et_bar_norm_rn', 'counts']]

print(test)

# think below is jsut for reference for the series, but not used in paper
columns = []
for i in range(1, 5):
  columns.append('term_%d' % i)
  mean_df[columns[-1]] = 1./(mean_df.lai*mean_df.uwue_norm*mean_df.g1**i)

print(mean_df.loc[:, columns])

mean_df['vpd_crit'] = et_min_vpd(mean_df, mean_df.lai)
mean_df['sigma_dot_uwue'] = mean_df.lai*mean_df.uwue_norm
columns = ['r_moist', 'c_a', 'gamma', 'lai',\
           'sigma_dot_uwue', 'vpd_crit']
print('\n vpd critical table')
print(mean_df.loc[:, columns])

def optimizer(vpds, *args):
  """this finds the otopmism vpd to get the most true hits"""
  _dff = args[0]
  count_true_less = np.array([float(_dff.loc[((_dff.vpd < vpd)\
                                    & (_dff.d_et < 0.0)), 'd_et'].count())\
                              for vpd in vpds])
  count_true_more = np.array([float(_dff.loc[((_dff.vpd > vpd)\
                                    & (_dff.d_et > 0.0)), 'd_et'].count())\
                              for vpd in vpds])
  output = (count_true_less+count_true_more)#/float(_dff.d_et.count())
  return output


def vpd_statistics(_df, ax, site=False):
  """computes fraction of time theory is correct, and mean det/dvpd as a
  function of theory"""
  columns = ['frac_correct_less', 'frac_incorrect_less',\
             'frac_correct_more', 'frac_incorrect_more',\
             'mean_less', 'mean_more', 'mean'] # ,\
             # 'optimized_vpd', 'true_vpd']
  # mean_df.loc[_df.pft.iloc[0], 'vpd_crit']
  median = _df.median()
  vpd_crit = et_min_vpd(_df.mean(), _df.mean().lai)#median.lai)
  less_df = _df.d_et.loc[(_df.vpd < vpd_crit)]
  try:
    less = [float(less_df.loc[less_df < 0.0].count())/float(less_df.count()),
            float(less_df.loc[less_df > 0.0].count())/float(less_df.count())]
  except ZeroDivisionError:
    less = [0.0, 0.0]

  more_df = _df.d_et.loc[(_df.vpd > vpd_crit)]
  try:
    more = [float(more_df.loc[more_df > 0.0].count())/float(more_df.count()),
            float(more_df.loc[more_df < 0.0].count())/float(more_df.count())]
  except ZeroDivisionError:
    more = [0.0, 0.0]
  # print('test', vpd_crit)
  # vpd_opt = scipy.optimize.golden(optimizer, args=(_df))# ,\
  #                                 # brack=(0.0, 5.0*vpd_crit),\
  # print(vpd_opt)
  vpds = np.linspace(0.0, _df.vpd.max())
  hits = optimizer(vpds, _df)
  if site:
    p = ax.plot(vpds, hits, label=_df.site.iloc[0])
    ax.set_title(_df.pft.iloc[0])
  else:
    p = ax.plot(vpds, hits, label=_df.pft.iloc[0])
  ax.plot([vpd_crit], [hits[np.absolute(vpds-vpd_crit).argmin()]],\
          marker='*', color=p[-1].get_color())
  #ax.plot([vpd_crit, vpd_crit], [hits.min(), hits.max()], label=_df.pft.iloc[0])

  data = np.hstack([less, more,\
                    [less_df.mean(), more_df.mean() , _df.d_et.mean()]])
  #, vpd_opt, vpd_crit]
  df_out = pd.DataFrame(data=[data], index=[_df.pft.iloc[0]], columns=columns)
  return df_out

fig = plt.figure()
ax = fig.add_subplot(111)
fraction_data = df.groupby('pft').apply(vpd_statistics, ax)
plt.legend(loc='best')
plt.savefig('./garb_fullcln_mean_all.png')
print(fraction_data)

def site_plot(pft):
  """just a wrapper for doing indiviual site plots"""
  fig = plt.figure()
  ax = fig.add_subplot(111)
  fraction_data = df.loc[(df.pft == pft), :]\
                    .groupby('site').apply(vpd_statistics, ax, site=True)
  plt.legend(loc='best')
  plt.savefig('./garb_%s_fullcln.png' % pft)
  print(fraction_data)
  return

# site_plot('DBF')
# site_plot('ENF')
# site_plot('CSH')
# site_plot('GRA')
# site_plot('CRO')


# df = pd.read_pickle('%s/changjie/full_pandas_fix_scaling.pkl'\
#                     % os.environ['DATA'])
# print(df.shape)
# df = df.groupby('site').apply(site_clean)
# print('site', df.shape)
# df = clean_df(df)
# print('cleaned', df.shape)
# df = clean_df(df, var='lai_gpp')
# print('lceaned gpp', df.shape)
# print(df.shape)


# meta = {}
# meta['folder_label'] = 'pft'
# meta['folder'] = 'hist_plots'
# meta['var'] = 'lai'
# meta['suff'] = 'site_clean'


# print('all removal')
# df = pd.read_pickle('%s/changjie/full_pandas_fix_scaling.pkl'\
#                     % os.environ['DATA'])
# print(df.shape)
# df = site_clean(df)
# print('outlier', df.shape)
# df = clean_df(df)
# print('cleaned', df.shape)
# df = clean_df(df, var='lai_gpp')
# print('lceaned gpp', df.shape)
# print(df.shape)


# meta = {}
# meta['folder_label'] = 'pft'
# meta['folder'] = 'hist_plots'
# meta['var'] = 'lai'
# meta['suff'] = 'full_clean'

# import codebase.plot_tools as plot_tools
# df.groupby('pft').apply(plot_tools.histogram, meta)

