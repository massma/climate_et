#! ~/edward/bin/python
"""
This script makes a map of PFT for the paper, fig1.
"""
import glob
import os
import pandas as pd
import numpy as np
import metcalcs as met
import seaborn as sns
import resource
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib as mpl
import codebase.penman_monteith as pm
import util
import importlib
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
mpl.rcParams.update(mpl.rcParamsDefault)

importlib.reload(util)

# grab sites actually used in analysis
df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA']).loc[:, 'site'].drop_duplicates()

# get metadata
meta = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                   'Site_list_(canopy_height).csv' % (os.environ['DATA']))
# subset by sites actually used
meta.index = meta.Site
meta = meta.loc[df.values]

def make_map(_df):
  """makes a map given a _df with Lat, Lon, and Cover"""
  pfts = _df.Cover_type.drop_duplicates()
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*3)
  axs = [fig.add_subplot(3, 1, i+1) for i in range(3)]
  xlims = [(-180., 180.), (-140., -50.), (-20., 40.)]
  ylims = [(-50., 70.), (15., 60.), (30., 60.)]
  for ax, xlim, ylim in zip(axs, xlims, ylims):
    print(ylim)
    print(xlim)
    m = Basemap(projection='merc', llcrnrlat=ylim[0], urcrnrlat=ylim[1],\
            llcrnrlon=xlim[0], urcrnrlon=xlim[1], lat_ts=np.mean(ylims),\
                resolution='c', ax=ax)#res was l
    m.drawcoastlines()
    for pft in pfts:
      print(pft)
      _ds = _df.loc[(_df.Cover_type==pft), ['Latitude', 'Longitude']]
      x, y = m(_ds.Longitude.values, _ds.Latitude.values)
      sc = ax.scatter(x, y, s=16, label=pft)
    # m.drawparallels(np.arange(-90.,120.,30.))
    # m.drawmeridians(np.arange(0.,420.,60.))
  plt.legend(loc='best')
  util.test_savefig('../doc/paper/figs/fig01.png')
  return

make_map(meta)


