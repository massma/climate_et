#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *
import scipy.optimize
import matplotlib.pyplot as plt


plt.close('all')
###### table 5 ####
def frequency(_df):
  """return fraction fos amples d_et < 0"""
  return _df.d_et[_df.d_et < 0.].count() / _df.d_et.count()


#### Table 1, pft tables
print('\n\n******TABLE 1*****')
print(mean_df.loc[:, ['g1', 'uwue', 'uwue_zhou', 'uwue_zhou_std']])


#### Table 2, pft tables
print('\n\n******TABLE 3*****')
print(mean_df.loc[:, ['r_moist', 'c_a', 'gamma', 'vpd_crit']])

std = df.groupby('pft').std()
mean_df['d_et_bar'] = d_et(mean_df)
mean_df['d_et_bar_std'] = np.absolute(mean_df['d_et_bar'])*std['vpd']
mean_df['d_et_d_rn_bar_std'] = mean_df['delta']/\
                               (mean_df['delta'] + mean_df['gamma'])\
                               *std['r_net']
mean_df['d_et_bar_norm'] = mean_df['d_et_bar_std']\
                           /std['et_obs']
mean_df['total_var'] = np.sqrt(mean_df.d_et_bar_std**2\
                               +mean_df.d_et_d_rn_bar_std**2)
mean_df['d_et_bar_norm_total'] = mean_df['d_et_bar_std']\
                           /mean_df['total_var']
mean_df['et_std'] = std['et_obs']
counts = df.groupby('pft').apply(frequency)

mean_df['counts'] = counts
# test = mean_df.loc[:, ['d_et', 'd_et_bar', 'd_et_bar_std',\
#                        'd_et_bar_norm', 'counts']]


print('\n\n******TABLE 4******')
print(mean_df.loc[:, ['d_et', 'd_et_bar', 'counts']])


print(mean_df.uwue)
print(mean_df.uwue_zhou)

# think below is jsut for reference for the series, but not used in paper
# columns = []
# for i in range(1, 5):
#   columns.append('term_%d' % i)
#   mean_df[columns[-1]] = 1./(mean_df.lai*mean_df.uwue_norm*mean_df.g1**i)

# print(mean_df.loc[:, columns])

# mean_df['vpd_crit'] = et_min_vpd(mean_df, mean_df.lai)
# mean_df['sigma_dot_uwue'] = mean_df.lai*mean_df.uwue_norm
# columns = ['r_moist', 'c_a', 'gamma', 'lai',\
#            'sigma_dot_uwue', 'vpd_crit']
# print('\n vpd critical table')
# print(mean_df.loc[:, columns])

# def optimizer(vpds, *args):
#   """this finds the otopmism vpd to get the most true hits"""
#   _dff = args[0]
#   count_true_less = np.array([float(_dff.loc[((_dff.vpd < vpd)\
#                                     & (_dff.d_et < 0.0)), 'd_et'].count())\
#                               for vpd in vpds])
#   count_true_more = np.array([float(_dff.loc[((_dff.vpd > vpd)\
#                                     & (_dff.d_et > 0.0)), 'd_et'].count())\
#                               for vpd in vpds])
#   output = (count_true_less+count_true_more)#/float(_dff.d_et.count())
#   return output


# def vpd_statistics(_df, ax, site=False):
#   """computes fraction of time theory is correct, and mean det/dvpd as a
#   function of theory"""
#   columns = ['frac_correct_less', 'frac_incorrect_less',\
#              'frac_correct_more', 'frac_incorrect_more',\
#              'total_frac_correct', 'mean_less', 'mean_more', 'mean'] # ,\
#              # 'optimized_vpd', 'true_vpd']
#   # mean_df.loc[_df.pft.iloc[0], 'vpd_crit']
#   # median = _df.median()
#   vpd_crit = et_min_vpd(_df.mean(), _df.mean().lai)#median.lai)
#   print(vpd_crit)
#   less_df = _df.d_et.loc[(_df.vpd < vpd_crit)]
#   try:
#     less = [float(less_df.loc[less_df < 0.0].count())/float(less_df.count()),
#             float(less_df.loc[less_df > 0.0].count())/float(less_df.count())]
#   except ZeroDivisionError:
#     less = [0.0, 0.0]

#   more_df = _df.d_et.loc[(_df.vpd > vpd_crit)]
#   try:
#     more = [float(more_df.loc[more_df > 0.0].count())/float(more_df.count()),
#             float(more_df.loc[more_df < 0.0].count())/float(more_df.count())]
#   except ZeroDivisionError:
#     more = [0.0, 0.0]
#   # print('test', vpd_crit)
#   # vpd_opt = scipy.optimize.golden(optimizer, args=(_df))# ,\
#   #                                 # brack=(0.0, 5.0*vpd_crit),\
#   # print(vpd_opt)
#   vpds = np.linspace(0.0, _df.vpd.max())
#   hits = optimizer(vpds, _df)
#   if site:
#     p = ax.plot(vpds, hits, label=_df.site.iloc[0])
#     ax.set_title(_df.pft.iloc[0])
#   else:
#     p = ax.plot(vpds, hits, label=_df.pft.iloc[0])
#   ax.plot([vpd_crit], [hits[np.absolute(vpds-vpd_crit).argmin()]],\
#           marker='*', color=p[-1].get_color())
#   #ax.plot([vpd_crit, vpd_crit], [hits.min(), hits.max()], label=_df.pft.iloc[0])

#   data = np.hstack([less, more,\
#                     [float(optimizer([vpd_crit], _df))/_df.d_et.count(),\
#                      less_df.mean(), more_df.mean() , _df.d_et.mean()]])
#   #, vpd_opt, vpd_crit]
#   df_out = pd.DataFrame(data=[data], index=[_df.pft.iloc[0]], columns=columns)
#   return df_out

# fig = plt.figure()
# ax = fig.add_subplot(111)
# fraction_data = df.groupby('pft').apply(vpd_statistics, ax)
# plt.legend(loc='best')
# ax.set_xlabel('VPD (Pa)')
# ax.set_ylabel('# observations theory gets right')
# plt.savefig('./for_agu.pdf')
# print(fraction_data)

# def site_plot(pft):
#   """just a wrapper for doing indiviual site plots"""
#   fig = plt.figure()
#   ax = fig.add_subplot(111)
#   fraction_data = df.loc[(df.pft == pft), :]\
#                     .groupby('site').apply(vpd_statistics, ax, site=True)
#   plt.legend(loc='best')
#   plt.savefig('./garb_%s_normalized.png' % pft)
#   print(fraction_data)
#   return

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

fig = plt.figure()
ax = fig.add_subplot(111)
for g1, uwue, label in zip(mean_df.g1, mean_df.uwue, mean_df.index):
  ax.scatter(g1, uwue, label=label)
ax.set_xlabel('g1')
ax.set_ylabel('uwue')
plt.legend(loc='best')
plt.savefig('%s/climate_et/g1_uwue_relationship.png' % os.environ['PLOTS'])
