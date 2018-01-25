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

def plot_box(_df, pft=False):
  """makes a box plot"""
  plt.close('all')
  fig = plt.figure()
  ax = fig.add_subplot(111)
  # exclude the last 5 percent b/c usually large
  ptiles = np.arange(0.0, 91.0, 5.0)
  positions = []
  labels = []
  data = []
  widths = []
  for ptile in ptiles:
    low = _df.swc.quantile(q=ptile/100.)
    high = _df.swc.quantile(q=ptile/100.+0.05)
    data.append(_df.uwue[(_df.swc > low) & (_df.swc <= high)])
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
  ax.set_ylabel('uWUE', fontsize=fontsize)
  ax.set_xlabel('SWC', fontsize=fontsize)
  if pft:
    ax.set_title('PFT: %s'\
                 % (_df.pft.iloc[0]),\
                 fontsize=fontsize+3)
  else:
    ax.set_title('Site: %s, PFT: %s'\
                 % (_df.site.iloc[0], _df.pft.iloc[0]),\
                 fontsize=fontsize+3)

  # ax.set_ylim(ylim)

  # ax.set_xlim(orig_xlim)
  plt.tight_layout()
  if pft:
    plt.savefig('%s/climate_et/swc_box/pft/%s_box.png'\
                % (os.environ['PLOTS'], _df.pft.iloc[0]))
  else:
    plt.savefig('%s/climate_et/swc_box/%s_%s_box.png'\
                % (os.environ['PLOTS'], _df.pft.iloc[0], _df.site.iloc[0]))
  return

#_df = df.iloc[:1000, :].copy()
#plot_box(_df)
#df.groupby('site').apply(plot_box)
df.groupby('pft').apply(plot_box, pft=True)

