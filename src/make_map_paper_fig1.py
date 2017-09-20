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
  ax = fig.add_subplot(111)
  m = Basemap(projection='npaeqd',boundinglat=10,lon_0=270,resolution='l')
  m.drawcoastlines()
  for pft in pfts:
    print(pft)
    _ds = _df.loc[(_df.Cover_type==pft), ['Latitude', 'Longitude']]
    x, y = m(_ds.Longitude.values, _ds.Latitude.values)
    sc = ax.scatter(x, y, s=50, label=pft)
  plt.legend(loc='best')
  util.test_savefig('../doc/paper/figs/fig01.png')
  return

make_map(meta)


