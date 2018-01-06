#! ~/edward/bin/python
"""
This script makes box plot figures for code
"""
from shared_functions import *
import seaborn as sns
import importlib

fontsize=18
linewidth=2.5
boxprops = {'linewidth' : 1.5}
medianprops = {'linewidth' : 2.5}

dashedlinewidth=1.4

bluelinewidth=3.0
# below is whether to plot taylor series approximation
PLOT_SERIES = False
#below is to make slight changes for the talk
PLOT_TALK = True

name_dict = {'CRO': 'Crops',\
             'DBF': 'Deciduous Forest',
             'ENF': 'Evergreen Forest',
             'GRA': 'Grass',
             'CSH': 'Shrub (closed)'}

paren_string = r'$\left(\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_s }{1.6 \; R\; \sigma \; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{VPD}}'\
               r'{2 (g_1 + \sqrt{VPD})^2}\right)\right)$'

def plot_box(_df):
  """makes a box plot"""
  mean_row = mean_df.loc[_df.pft.iloc[0], :]
  plt.close('all')
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ptiles = np.arange(5.0, 91.0, 5.0)
  positions = []
  labels = []
  data = []
  widths = []
  for ptile in ptiles:
    low = _df.vpd.quantile(q=ptile/100.)
    high = _df.vpd.quantile(q=ptile/100.+0.05)
    data.append(_df.sign[(_df.vpd > low) & (_df.vpd <= high)])
    positions.append((low+high)*0.5/1.0e3)
    widths.append((high-low)*0.8/1.0e3)
    labels.append('%4.0f' % positions[-1])
  ax.boxplot(x=data, positions=positions,\
             sym=(''), widths=widths, manage_xticks=False,\
             whis=[5, 95], boxprops=boxprops, medianprops=medianprops)#sym=('')
  print('make hist for %s' % _df.pft.iloc[0])
  orig_xlim = ax.get_xlim()
  #ax.plot(orig_xlim, [d_calc.CP/287.0, d_calc.CP/287.0], 'm--', linewidth=1.0)
  ylim = ax.get_ylim()#[-4.0, 4.0]
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  uwue = mean_row.uwue
  #\ uwue=mean_row.uwue_zhou),\
  ax.plot(vpd/1.e3, d_calc.sign(mean_row, vpd=vpd),
          'b-', linewidth=bluelinewidth)
  # ax.plot(vpd/1.e3, d_calc.sign(mean_row, vpd=vpd, uwue=(mean_row.uwue_zhou\
  #                                             +mean_row.uwue_zhou_std)),\
  #         'b-', linewidth=1.0)
  # ax.plot(vpd/1.e3, d_calc.sign(mean_row, vpd=vpd, uwue=(mean_row.uwue_zhou\
  #                                               -mean_row.uwue_zhou_std)),\
  #         'b-', linewidth=1.0)
  #_df.groupby(pd.cut(_df.vpd, bins=20)).apply(add_box, ax=ax)
  #_df.boxplot(column='sign', by=pd.cut(_df.vpd, bins=20), ax=ax)
  if (_df.pft.iloc[0] == 'ENF') | (_df.pft.iloc[0] == 'DBF')\
     | (_df.pft.iloc[0] == 'CSH'):
    ax.set_ylabel(paren_string, fontsize=fontsize)
  if (_df.pft.iloc[0] == 'CSH') | (_df.pft.iloc[0] == 'GRA'):
    ax.set_xlabel('VPD (kPa)', fontsize=fontsize)
  ax.set_title(name_dict[_df.pft.iloc[0]], fontsize=fontsize+3)
  ax.plot(orig_xlim, [0., 0.], 'k--', linewidth=dashedlinewidth)

  ax.set_ylim(ylim)

  ax.set_xlim(orig_xlim)
  # if _df.pft.iloc[0] == 'ENF':
  #   ax.text(1.95, 1.2, '*', horizontalalignment='center',\
  #           verticalalignment='center', fontdict={'fontsize' : 40})
  plt.tight_layout()
  plt.savefig('../../doc/shared_figs/%s_box.pdf' % _df.pft.iloc[0])
  return

#_df = df.iloc[:1000, :].copy()
#plot_box(_df)
df.groupby('pft').apply(plot_box)

os.system('pdfjam %s/DBF_box.pdf %s/ENF_box.pdf %s/CSH_box.pdf '\
          '--nup 1x3 --no-landscape --outfile %s/first-column.pdf'\
          % tuple(['../../doc/shared_figs']*4))
os.system('pdfcrop --margins 10 %s %s'\
          % tuple(['../../doc/shared_figs/first-column.pdf']*2))
os.system('pdfjam %s/CRO_box.pdf %s/GRA_box.pdf '\
          '--nup 1x3 --no-landscape --outfile %s/second-column.pdf'\
          % tuple(['../../doc/shared_figs']*3))
os.system('pdfcrop --margins 10 %s %s'\
          % tuple(['../../doc/shared_figs/second-column.pdf']*2))
os.system('pdfjam %s %s '\
          '--nup 2x1 --no-landscape --outfile %s/agu_fig.pdf'\
          % ('../../doc/shared_figs/first-column.pdf',\
             '../../doc/shared_figs/second-column.pdf',\
             '../../doc/paper/'))
os.system('pdfcrop --margins 10 %s %s'\
          % tuple(['../../doc/paper/agu_fig.pdf']*2))
