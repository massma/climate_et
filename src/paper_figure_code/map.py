#! ~/edward/bin/python
"""
This script makes fig 01 for the paper
"""
from shared_functions import *

symbols = ['o', 'v', '>', 's', '<', 'D', 'X', 'P', '8']
plt.close('all')
def make_map(_df):
  """makes a map given a _df with Lat, Lon, and Cover"""
  pfts = _df.Cover_type.drop_duplicates()
  fig = plt.figure()
  #fig.set_figheight(fig.get_figheight())
  ax = fig.add_subplot(111)
  m = Basemap(projection='merc', llcrnrlat=-40, urcrnrlat=70,\
            llcrnrlon=-180, urcrnrlon=180, lat_ts=0.0,\
                resolution='l', ax=ax)#res was l
  m.fillcontinents(alpha=1.0)
  size = 10

  for i,pft in enumerate(pfts):
    print(pft)
    _ds = _df.loc[(_df.Cover_type == pft), ['Latitude', 'Longitude']]
    x, y = m(_ds.Longitude.values, _ds.Latitude.values)
    ax.scatter(x, y, s=12, label=pft, alpha=0.7,\
               marker=symbols[i], zorder=99999)
  ax.set_title('FLUXNET2015 Sites with $\geq4$ Years Data')


  plt.legend(loc='best')
  plt.tight_layout()
  util.test_savefig('../../doc/paper/map.pdf')
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
