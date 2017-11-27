#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""

from shared_functions import *

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

