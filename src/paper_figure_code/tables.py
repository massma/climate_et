#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
from shared_functions import *
import scipy.optimize
import fluxnet_pycite.fluxnet_pycite as fluxcite

importlib.reload(fluxcite)

plt.close('all')

##### fluxnet citation table"""
fluxcite.generate_table(list(sites_used), "../../doc/paper/flux_sites.tex",
                        "../../doc/paper/flux_sites.bib")


###### table 5 ####
def frequency(_df):
  """return fraction fos amples d_et < 0"""
  return _df.d_et[_df.d_et < 0.].count() / _df.d_et.count()

counts = df.groupby('pft').apply(frequency)
mean_df['counts'] = counts

pft_order = mean_df.sort_values('counts').index

#### Table 1, pft tables
print('\n\n******TABLE 1*****')
mean_df['g1-dot-uwue'] = mean_df.g1*mean_df.uwue
print(mean_df.loc[:, ['g1', 'uwue', 'uwue_zhou', 'uwue_zhou_std']])#,\
                      #'g1-dot-uwue']])
def write_table_1(_df):
  """write output of table 1"""
  fh = open("../../doc/paper/pft_params.tex", "w")
  for pft in pft_order:
    row = _df.loc[pft, :]
    if (row.uwue_zhou == 1.0) & (row.uwue_zhou_std == 1.0):
      fh.write("%s & %s & %.2f & %.2f & N/A \\\\ \n"
               % (pft, name_dict[pft][:-5], np.round(row.g1, 2),
                  np.round(row.uwue, 2)))
    else:
      fh.write("%s & %s & %.2f & %.2f & %.2f $\pm$ %.2f \\\\ \n"
               % (pft, name_dict[pft][:-5], np.round(row.g1, 2),
                  np.round(row.uwue, 2), np.round(row.uwue_zhou, 2),
                  np.round(row.uwue_zhou_std, 2)))

  fh.close()
  return

write_table_1(mean_df)

#### Table 2, pft tables
print('\n\n******TABLE 3*****')
def write_table_2(_df):
  """write output of table 1"""
  fh = open("../../doc/paper/vpd_crit.tex", "w")
  for pft in pft_order:
    row = _df.loc[pft, :]
    fh.write("%s & %.1f & %.1f & %.1f & \\textbf{%.1f}  \\\\ \n"
             % (pft,  np.round(row.r_moist, 1),
                np.round(row.c_a, 1), np.round(row.gamma, 1),
                np.round(row.vpd_crit, 1)))

  fh.close()
  return

write_table_2(mean_df)

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


print('\n\n******TABLE 4******')
print(mean_df.loc[:, ['d_et', 'd_et_bar', 'counts']])
def write_table_3(_df):
  """write output of table 1"""
  fh = open("../../doc/paper/stats.tex", "w")
  for pft in pft_order:
    row = _df.loc[pft, :]
    fh.write("%s & %.3f & %.3f & %.3f \\\\ \n"
             % (pft,  np.round(row.d_et, 3),
                np.round(row.d_et_bar, 3), np.round(row.counts, 3)))

  fh.close()
  return

write_table_3(mean_df)


fig = plt.figure()
ax = fig.add_subplot(111)
for g1, uwue, label in zip(mean_df.g1, mean_df.uwue, mean_df.index):
  ax.scatter(g1, uwue, label=label)
ax.set_xlabel('g1')
ax.set_ylabel('uwue')
plt.legend(loc='best')
plt.savefig('%s/climate_et/g1_uwue_relationship.png' % os.environ['PLOTS'])
