#! ~/edward/bin/python
"""
This script makes all figs for the paper
"""
import os
import importlib
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import codebase.plot_tools as plot_tools
import util
import metcalcs as met
import codebase.penman_monteith as pm
import codebase.calc_tools as calc

mpl.use('Pdf')
mpl.rcParams.update(mpl.rcParamsDefault)
importlib.reload(util)
importlib.reload(plot_tools)
# grab sites actually used in analysis
df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA']).loc[:, 'site'].drop_duplicates()

# get metadata
site_list = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                   'Site_list_(canopy_height).csv' % (os.environ['DATA']))
# subset by sites actually used
site_list.index = site_list.Site
site_list = site_list.loc[df.values]

# below is the label 
paren_string = r'(Term 2 - Term 3) $\left(\frac{ c_p}{R_{air}} '\
               r'- \frac{\gamma c_s }{LAI \; 1.6 \; R\; uWUE  }'\
               r'\left( \frac{2 g_1 + \sqrt{D}}'\
               r'{2 (g_1 + \sqrt{D})^2}\right)\right)$'
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
  util.test_savefig('../doc/paper/fig01.pdf')
  return

### FIGURE 1 ###
#make a map fig 1
make_map(site_list)

### FIGURE 2 ###
#for final npaper should make this .pdf, and probably
df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])
df['g_a'] = 1./df['r_a']
df['d_et_leaf'] = df['scaling']*df['vpd_leaf']
df['d_et_atm'] = df['scaling']*df['vpd_atm']
df['uwue_norm'] = df.uwue/pm.LV

meta = {}
meta['var'] = 'lai'
meta['folder'] = ''
meta['folder_label'] = ''
plot_tools.histogram(df, meta)
os.system('cp %s/climate_et/histogram/full_lai.png '\
          '../doc/paper/fig02.png' % (os.environ['PLOTS']))

# def make_quant_map(_df, meta):
#   """makes a map given a _df with Lat, Lon, and Cover"""
#   fig = plt.figure()
#   fig.set_figheight(fig.get_figheight()*4)
#   axs = [fig.add_subplot(4, 1, i+1) for i in range(4)]
#   xlims = [(-140., -50.), (-20., 40.), (10., 40), (110., 155.)]
#   ylims = [(15., 60.), (30., 70.), (-20., -10.), (-40., -10.)]
#   sizes = [24, 12, 40, 40]
#   for ax, xlim, ylim, size in zip(axs, xlims, ylims, sizes):
#     print(ylim)
#     print(xlim)
#     m = Basemap(projection='merc', llcrnrlat=ylim[0], urcrnrlat=ylim[1],\
#             llcrnrlon=xlim[0], urcrnrlon=xlim[1], lat_ts=np.mean(ylims),\
#                 resolution='l', ax=ax)#res was l
#     m.drawcoastlines()
#     vmax = 0.07
#     x, y = m(_df.lon.values, _df.lat.values)
#     color = ax.scatter(x, y, c=_df[meta['var']], alpha=0.5,\
#                     s=size, cmap='RdBu',\
#                     vmin=-vmax, vmax=vmax)
#     # m.drawparallels(np.arange(-90.,120.,30.))
#     # m.drawmeridians(np.arange(0.,420.,60.))
#   plt.colorbar(color)
#   plt.tight_layout()
#   util.test_savefig('../doc/paper/fig03a.pdf')
#   return

# #create figure 2
# meta = {}
# meta['var'] = 'd_et'
# mean_site = df.groupby('site').mean()
# mean_site['lat'] = site_list.loc[mean_site.index, 'Latitude']
# mean_site['lon'] = site_list.loc[mean_site.index, 'Longitude']
# make_quant_map(mean_site, meta)

# plt.figure()
# x = np.linspace(0., 1., 10)
# plt.plot(x, x, 'k,')
# plt.scatter(x, x+0.1, s=1., c='k', edgecolor='')
# plt.savefig('%s/temp/test_pixel.png' % os.environ['PLOTS'])

### FIGURE 3 ###
# joint distribution of lai and swc
meta = {}
meta['ylim'] = (0., 2.)
meta['xlim'] = (0., 4000.)
meta['plot_type'] = '' #'simple'
meta['x_var'] = 'vpd'
meta['y_var'] = 'lai'
plot_tools.scatter_wrapper(df, meta)
os.system('cp %s/climate_et/scatters/vpd_lai.png ../doc/paper/fig03.png'\
          % (os.environ['PLOTS']))

### FIGURE 4 ###
# scaling term as a function of t and g_a
for key in ['gamma', 'delta', 'rho_a', 'g_a']:
  print('%s var: %f, mean: %f, cf: %f'\
        % (key, df[key].std(), df[key].mean(), df[key].std()/df[key].mean()))
df['rho_like'] = df['p_a']/(273.15 + df['t_a'])

def scaling(_df, t_a, g_a):
  """
  calculates the scaling term (Term 1 in paper) given _df, g_a (contant)
  and t_a
  """
  _atmos = {'t_a' : t_a, 'e_s' : met.vapor_pres(t_a)*pm.VP_FACTOR}
  delta = pm.delta(_atmos)
  return g_a*_df['rho_like'].mean()/(pm.delta(_atmos) + df.gamma.mean())

def plot_scaling(_df, ax, savefig=False):
  """makes idealized plots of scaling as a function of g_a and T"""
  t_a = np.linspace(_df.t_a.quantile(q=0.05), _df.t_a.quantile(q=0.95))
  ax.plot(t_a, np.ones(t_a.shape)*_df.scaling.mean(), 'k--',\
          linewidth=0.5, label='Term 1 mean')
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    g_a = _df.g_a.quantile(q=percentile/100.)
    scale = scaling(_df, t_a, g_a)
    ax.plot(t_a, scale, label='$g_a$ = %5.3f (%dth percentile)'\
            % (g_a, int(percentile)))
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Term 1 ($\frac{g_a \; P}{T(\Delta + \gamma)}$)')
  ax.set_title('PFT = %s' % _df.pft.iloc[0])
  plt.legend(loc='best', fontsize=8)
  if savefig:
    plt.savefig('../doc/paper/fig04.pdf')
  return

def scaled_mean(_df):
  """just scales g_a by mean"""
  _df['scaled_g_a'] = _df.g_a/df.g_a.mean()
  return _df

def plot_scaled_scaling(_df, ax):
  """makes idealized plots of scaling as a function of g_a and T"""
  t_a = np.linspace(_df.t_a.quantile(q=0.05), _df.t_a.quantile(q=0.95))
  _df = _df.groupby('pft').apply(scaled_mean)
  for percentile in [5., 25., 50., 75., 95.]:
    g_a = _df.scaled_g_a.quantile(q=percentile/100.)
    scale = scaling(_df, t_a, g_a)
    ax.plot(t_a, scale,\
            label=r'$\frac{g_a}{\overline{g_a}}$ = %5.3f (%dth percentile)'\
            % (g_a, int(percentile)))
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Normalized Term 1 '\
                r'($\frac{g_a \; P}{\overline{g_a} \; T(\Delta + \gamma)}$)')
  plt.legend(loc='best', fontsize=8)
  return

def plot_mean_pft(_df, ax, savefig=False):
  """makes idealized plots of scaling as a function of pft"""

  # pfts = ['GRA', 'ENF', 'CRO', 'DBF', 'CSH']
  pfts = ['DBF', 'ENF', 'CSH', 'CRO', 'GRA']
  for pft in pfts:

    idx = (_df.pft == pft)
    t_a = np.linspace(_df.t_a.loc[idx].quantile(q=0.05),\
                      _df.t_a.loc[idx].quantile(q=0.95))
    g_a = _df.g_a.loc[idx].mean()
    scale = scaling(_df, t_a, g_a)
    ax.plot(t_a, scale, label='PFT = %s, $\overline{g_a}$ = %5.3f'\
            % (pft, g_a))
    ptiles = np.array([_df.t_a.loc[idx].quantile(q=_p/100.)\
                       for _p in [25., 50., 75.]])
    ax.plot(ptiles, scaling(_df, ptiles, g_a), 'k*')
  ax.set_xlabel('T (C)')
  ax.set_ylabel(r'Term 1 ($\frac{g_a \; P}{T(\Delta + \gamma)}$)')
  plt.legend(loc='best', fontsize=8)
  return

def scaling_wrapper(df):
  """wrapper that groups by pft and does scaling plot"""

  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*2)
  ax = fig.add_subplot(2, 1, 1)
  plot_scaled_scaling(df, ax)
  ax = fig.add_subplot(2, 1, 2)
  plot_mean_pft(df, ax)
  plt.tight_layout()
  plt.savefig('../doc/paper/fig04.pdf')
  return


# plot_scaling(df)
# os.system('convert -append %s/climate_et/scaling/*.pdf '\
#           '../doc/paper/fig04.pdf'\
#           % (os.environ['PLOTS']))
scaling_wrapper(df)




##### Figure 5 #####
def term_2(_df, lai, vpd):
  """calculates term 2"""
  atmos = {'gamma' : _df.gamma.mean(), 'c_a' : _df.c_a.mean(),\
           'vpd' : vpd}
  if _df.uwue.std() > 1.e-8:
    print('error, uWUE is variable: %f!!!!' % _df.uwue.std())
  elif _df.g1.std() > 1.e-8:
    print('error, g1 is variabile: %f!!!!!' % _df.g1.std())
  canopy = {'uwue' : _df.uwue.mean(), 'g1' : _df.g1.mean()}
  return pm.CP/_df.r_moist.mean() +calc.leaf_vpd(atmos, canopy, lai)

def plot_leaf_vpd(_df, ax, savefig=False):
  """makes idealized plots of plant term as a function of lai and vpd"""
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    lai = _df.lai.quantile(q=percentile/100.)
    ax.plot(vpd, term_2(_df, lai, vpd), label='$LAI$ = %5.2f (%dth percentile)'\
            % (lai, int(percentile)))
  ax.set_xlabel('VPD (Pa)')
  ax.set_ylabel(paren_string)
  ax.set_title('PFT = %s, uWUE = %5.2f, g1 = %5.2f'\
               % (_df.pft.iloc[0], _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
  ax.plot(ax.get_xlim(), [0., 0.], 'k-', linewidth=0.2)
  plt.legend(loc='best', fontsize=8)
  if savefig:
    plt.savefig('../doc/paper/fig05.pdf')
  return


def plot_leaf_lai(_df, ax, savefig=False):
  """makes idealized plots of plant term as a function of lai and vpd"""
  lai = np.linspace(_df.lai.quantile(q=0.05), _df.lai.quantile(q=0.95))
  for percentile in [5., 25., 50., 75., 95.][::-1]:
    vpd = _df.vpd.quantile(q=percentile/100.)
    ax.plot(lai, term_2(_df, lai, vpd),\
            label='$VPD$ = %5.0f Pa (%dth percentile)'\
            % (vpd, int(percentile)))
  ax.set_xlabel('LAI')
  ax.set_ylabel(paren_string)
  ax.set_title('PFT = %s, uWUE = %5.2f, g1 = %5.2f'\
               % (_df.pft.iloc[0], _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
  ax.plot(ax.get_xlim(), [0., 0.], 'k-', linewidth=0.2)
  plt.legend(loc='best', fontsize=8)
  if savefig:
    plt.savefig('../doc/paper/fig05.pdf')
  return

def leaf_wrapper(df):
  """wraps df"""
  pfts = ['DBF', 'ENF', 'CSH', 'CRO', 'GRA']
  nplots = len(pfts)
  fig = plt.figure()
  fig.set_figheight(fig.get_figheight()*nplots)
  fig.set_figwidth(fig.get_figwidth()*2)
  print('\n')
  for i, pft in enumerate(pfts):
    _df = df.loc[(df.pft == pft), :]
    print('for %s, g1: %f, uwue: %f, vpd: %f, lai: %f'\
          % (pft, _df.g1.mean(), _df.uwue.mean()/pm.LV,\
             _df.vpd.mean(), _df.lai.mean()))
    ax = fig.add_subplot(nplots, 2, i*2+1)
    ax2 = fig.add_subplot(nplots, 2, i*2+2)
    plot_leaf_vpd(_df, ax)
    plot_leaf_lai(_df, ax2)
  plt.tight_layout()
  plt.savefig('../doc/paper/fig05.pdf')
  return

plt.close('all')
leaf_wrapper(df)
grouped = df.groupby('pft')
print('mean lai: %5.2f' % grouped.lai.mean().mean())
print('mean vpd: %5.2f' % grouped.vpd.mean().mean())

for key in ['lai', 'vpd', 'g1', 'uwue_norm']:
  print('cv %s: %5.2f' % (key,\
                          grouped[key].mean().std()/grouped[key].mean().mean()))

def get_pft(_df):
  return _df['pft'].iloc[0]

names = grouped.apply(get_pft)
plt.figure()
plt.plot(grouped.g1.mean(), grouped.uwue.mean(), 'k*')
for name, x, y in zip(names, grouped.g1.mean(), grouped.uwue.mean()):
  plt.annotate(name, xy=(x, y))
plt.savefig('%s/temp/garb.png' % os.environ['PLOTS'])

#now for figure 6 split product into two

### FIGURE 6 ###
# look at what controls variability between pfts
def first_half(_df, lai):
  """calcs first half of term 3"""
  return -_df.gamma.mean()*df.c_a.mean()/\
(lai*1.6*pm.R_STAR*_df.uwue_norm.mean())

def second_half(_df, vpd):
  """calcs second half of term 3"""
  _g1 = _df.g1.iloc[0]
  return -(2.*_g1 + np.sqrt(vpd))/\
    (2.*(_g1 + np.sqrt(vpd))**2)

def pft_leaf(_df, axs):
  """takes df and plots both halves of product in term 2"""
  lai = np.linspace(_df.lai.quantile(q=0.05), _df.lai.quantile(q=0.95))
  vpd = _df.vpd.mean()
  axs[0].plot(lai, term_2(_df, lai, vpd),\
              label=r'%s, $\overbar{VPD}$ = %4.0f Pa'\
              % (_df.pft.iloc[0], vpd))
  axs[1].plot(lai, first_half(_df, lai), label=str(_df.pft.iloc[0]))
  ptiles = np.array([_df.lai.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(ptiles, first_half(_df, ptiles), 'k*')
  # now second half
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  lai = _df.lai.mean()
  axs[2].plot(vpd, term_2(_df, lai, vpd),\
              label=r'%s, $\overbar{LAI}$ = %4.2f'\
              % (_df.pft.iloc[0], lai))
  axs[3].plot(vpd, second_half(_df, vpd), label=str(_df.pft.iloc[0]))
  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[3].plot(ptiles, second_half(_df, ptiles), 'k*')
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*2)
fig.set_figwidth(fig.get_figwidth()*2)
axs = [fig.add_subplot(2, 2, i+1) for i in range(4)]
df.groupby('pft').apply(pft_leaf, axs)
axs[0].set_xlabel('LAI')
axs[1].set_xlabel('LAI')
axs[0].set_ylabel(paren_string)
axs[1].set_ylabel(r'-$\frac{\gamma c_s }{LAI \; 1.6 \; R\;  uWUE }$')
axs[2].set_xlabel('VPD (Pa)')
axs[3].set_xlabel('VPD (Pa)')
axs[2].set_ylabel(paren_string)
axs[3].set_ylabel(r'-$\frac{2 g_1 + \sqrt{D}}{2 (g_1 + \sqrt{D})^2}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../doc/paper/fig06.pdf')
