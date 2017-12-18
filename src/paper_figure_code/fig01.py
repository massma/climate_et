#! ~/edward/bin/python
"""
This script makes fig 01 for the paper
"""
from shared_functions import *

def make_map(_df):
  """makes a map given a _df with Lat, Lon, and Cover"""
  pfts = _df.Cover_type.drop_duplicates()
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*4)
  axs = [fig.add_subplot(4, 1, i+1) for i in range(4)]
  xlims = [(-140., -50.), (-20., 40.), (10., 40), (110., 155.)]
  ylims = [(15., 60.), (30., 70.), (-20., -10.), (-40., -10.)]
  sizes = [24, 12, 40, 40]
  for ax, xlim, ylim, size in zip(axs, xlims, ylims, sizes):
    print(ylim)
    print(xlim)
    m = Basemap(projection='merc', llcrnrlat=ylim[0], urcrnrlat=ylim[1],\
            llcrnrlon=xlim[0], urcrnrlon=xlim[1], lat_ts=np.mean(ylims),\
                resolution='l', ax=ax)#res was l
    m.drawcoastlines()
    for pft in pfts:
      print(pft)

      _ds = _df.loc[(_df.Cover_type == pft), ['Latitude', 'Longitude']]
      x, y = m(_ds.Longitude.values, _ds.Latitude.values)
      ax.scatter(x, y, s=size, label=pft, alpha=1.)
    # m.drawparallels(np.arange(-90.,120.,30.))
    # m.drawmeridians(np.arange(0.,420.,60.))
  plt.legend(loc='best')
  plt.tight_layout()
  util.test_savefig('../../doc/paper/fig01_garb.pdf')
  return

### FIGURE 1 ###
#make a map fig 1
# get metadata
sites_used = df.loc[:, 'site'].drop_duplicates()

site_list = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                   'Site_list_(canopy_height).csv' % (os.environ['DATA']))
site_list.index = site_list.Site
site_list = site_list.loc[sites_used.values]

make_map(site_list)
