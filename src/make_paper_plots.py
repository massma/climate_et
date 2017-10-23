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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import codebase.plot_tools as plot_tools
import util
import metcalcs as met
import codebase.penman_monteith as pm
import codebase.calc_tools as calc
from sympy import Symbol, sqrt, series, latex, limit


#try some analytis with sympy
g1 = Symbol('g_1')
x = Symbol('\frac{\sqrt{D}}{g_1}')
func = g1*(2 + x)/(2*g1**2*(1 + x)**2)
print(latex(series(func, x, x0=0., dir='+', n=4)))
# x = Symbol('\sqrt{D}')
# func  = (2*g1 + x)/(2*(g1 + x)**2)
# print(latex(series(func, x, x0=0., dir='+')))

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
df['r_net'] = df['r_n'] - df['g_flux']
jacobians = {'vpd' : calc.d_et,\
             'lai' : calc.d_et_d_lai,\
             'seasonal_lai' : calc.d_et_d_lai,\
             'residual_lai' : calc.d_et_d_lai,\
             'g_a' : calc.d_et_d_g_a,\
             'delta' : calc.d_et_d_delta,\
             'r_net' : calc.d_et_d_r_net}

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
#   plt.tight_layout()n
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




# ##### Figure 5 #####
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


def term_2_approx(_df, lai, vpd, order=4):
  """calculates term 2"""
  atmos = {'gamma' : _df.gamma.mean(), 'c_a' : _df.c_a.mean(),\
           'vpd' : vpd}
  if _df.uwue.std() > 1.e-8:
    print('error, uWUE is variable: %f!!!!' % _df.uwue.std())
  elif _df.g1.std() > 1.e-8:
    print('error, g1 is variabile: %f!!!!!' % _df.g1.std())
  canopy = {'uwue' : _df.uwue.mean(), 'g1' : _df.g1.mean()}
  if order == 4:
    return pm.CP/_df.r_moist.mean() \
      -atmos['gamma']*atmos['c_a']*pm.LV/\
      (lai*1.6*pm.R_STAR*canopy['uwue'])\
      *(1./canopy['g1'] - 3.*np.sqrt(atmos['vpd'])/(2.*canopy['g1']**2)
        + 2.*atmos['vpd']/canopy['g1']**3\
        - 5.*np.sqrt(atmos['vpd'])**3/(2.*canopy['g1']**4))
  elif order == 2:
    return pm.CP/_df.r_moist.mean() \
      -atmos['gamma']*atmos['c_a']*pm.LV/\
      (lai*1.6*pm.R_STAR*canopy['uwue'])\
      *(1./canopy['g1'] - 3.*np.sqrt(atmos['vpd'])/(2.*canopy['g1']**2))
  else:
    print('error uncoded order number %d' % order)
    return

# def plot_leaf_vpd(_df, ax, savefig=False):
#   """makes idealized plots of plant term as a function of lai and vpd"""
#   vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
#   for percentile in [5., 25., 50., 75., 95.][::-1]:
#     lai = _df.lai.quantile(q=percentile/100.)
#     ax.plot(vpd, term_2(_df, lai, vpd), label='$LAI$ = %5.2f (%dth percentile)'\
#             % (lai, int(percentile)))
#   ax.set_xlabel('VPD (Pa)')
#   ax.set_ylabel(paren_string)
#   ax.set_title('PFT = %s, uWUE = %5.2f, g1 = %5.2f'\
#                % (_df.pft.iloc[0], _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
#   ax.plot(ax.get_xlim(), [0., 0.], 'k-', linewidth=0.2)
#   plt.legend(loc='best', fontsize=8)
#   if savefig:
#     plt.savefig('../doc/paper/fig05.pdf')
#   return


# def plot_leaf_lai(_df, ax, savefig=False):
#   """makes idealized plots of plant term as a function of lai and vpd"""
#   lai = np.linspace(_df.lai.quantile(q=0.05), _df.lai.quantile(q=0.95))
#   for percentile in [5., 25., 50., 75., 95.][::-1]:
#     vpd = _df.vpd.quantile(q=percentile/100.)
#     ax.plot(lai, term_2(_df, lai, vpd),\
#             label='$VPD$ = %5.0f Pa (%dth percentile)'\
#             % (vpd, int(percentile)))
#   ax.set_xlabel('LAI')
#   ax.set_ylabel(paren_string)
#   ax.set_title('PFT = %s, uWUE = %5.2f, g1 = %5.2f'\
#                % (_df.pft.iloc[0], _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
#   ax.plot(ax.get_xlim(), [0., 0.], 'k-', linewidth=0.2)
#   plt.legend(loc='best', fontsize=8)
#   if savefig:
#     plt.savefig('../doc/paper/fig05.pdf')
#   return

# def leaf_wrapper(df):
#   """wraps df"""
#   pfts = ['DBF', 'ENF', 'CSH', 'CRO', 'GRA']
#   nplots = len(pfts)
#   fig = plt.figure()
#   fig.set_figheight(fig.get_figheight()*nplots)
#   fig.set_figwidth(fig.get_figwidth()*2)
#   print('\n')
#   for i, pft in enumerate(pfts):
#     _df = df.loc[(df.pft == pft), :]
#     print('for %s, g1: %f, uwue: %f, vpd: %f, lai: %f'\
#           % (pft, _df.g1.mean(), _df.uwue.mean()/pm.LV,\
#              _df.vpd.mean(), _df.lai.mean()))
#     ax = fig.add_subplot(nplots, 2, i*2+1)
#     ax2 = fig.add_subplot(nplots, 2, i*2+2)
#     plot_leaf_vpd(_df, ax)
#     plot_leaf_lai(_df, ax2)
#   plt.tight_layout()
#   plt.savefig('../doc/paper/fig05.pdf')
#   return

# plt.close('all')
# leaf_wrapper(df)

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

def et_min_vpd(_df, lai):
  """calculates theoretical max vpd as functoin of -df and lai"""
  c3 = pm.CP/_df.r_moist
  c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
  c2 = _df.g1
  sqrt_vpd = (c1 + np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3)
  try:
    sqrt_vpd[sqrt_vpd < 0.] = np.nan
  except TypeError:
    if sqrt_vpd < 0.:
      sqrt_vpd = np.nan
  return sqrt_vpd**2

# def et_min_vpd1(_df, lai):
#   """calculates theoretical max vpd as functoin of -df and lai"""
#   """note below is only valid for negative x, which we don't have"""
#   c3 = pm.CP/_df.r_moist
#   c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
#   c2 = _df.g1
#   return ((c1 - np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3))**2


I = 5
def pft_leaf(_df, axs):
  """takes df and plots both halves of product in term 2"""
  global I
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  lai = _df.lai.mean()
  p = axs[0].plot(vpd, term_2(_df, lai, vpd),\
                  label="%s: $\overline{LAI}$=%4.2f, \n uWUE=%4.2f, g1=%4.1f"\
                  % (_df.pft.iloc[0],\
                     lai, _df.uwue_norm.iloc[0],  _df.g1.iloc[0]))
  axs[0].plot(vpd, term_2_approx(_df, lai, vpd, order=4), linestyle='dashed',\
              color=p[0].get_color())
  #print('axlim: ',axs[0].get_ylim())
  # axs[3].plot(vpd, second_half(_df, vpd),\
  #             label='%s, g1 = %4.1f' % (_df.pft.iloc[0], _df.g1.iloc[0]))
  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(vpd, np.ones(vpd.shape)*I)
  axs[1].plot(ptiles, np.ones(ptiles.shape)*I, 'k*')
  axs[-1].plot(np.ones(vpd.shape)*I, vpd)
  axs[-1].plot(np.ones(ptiles.shape)*I, ptiles, 'k*')

  # # axs[2].plot(ptiles, term_2(_df, lai, ptiles), 'k*')
  # axs[3].plot(ptiles, second_half(_df, ptiles), 'k*')

  lai = np.linspace(_df.lai.quantile(q=0.05), _df.lai.quantile(q=0.95))
  vpd = _df.vpd.mean()
  _mean_df = _df.mean()
  et_min = et_min_vpd(_mean_df, lai)
  # et_min[et_min > _df.vpd.quantile(0.99)] = np.nan
  # et_min[et_min > 4001.] = np.nan
  axs[2].plot(lai, et_min,\
              label=r"PFT = %s, uWUE=%4.2f, g1=%4.1f"\
              % (_df.pft.iloc[0],\
                 _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
  ptiles = np.array([_df.lai.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[3].plot(lai, np.ones(lai.shape)*I)
  axs[3].plot(ptiles, np.ones(ptiles.shape)*I, 'k*')
  # axs[1].plot(lai, first_half(_df, lai),\
  #             label='%s, uWUE = %4.2f'\
  #             % (_df.pft.iloc[0], _df.uwue_norm.iloc[0]))
  # # axs[0].plot(ptiles, term_2(_df, ptiles, vpd), 'k*')
  # axs[1].plot(ptiles, first_half(_df, ptiles), 'k*')
  # now second half
  I -= 1
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*2.5)
#fig.set_figwidth(fig.get_figwidth()*2)
#axs = [fig.add_subplot(2, 2, i+1) for i in range(4)]
_axs = [fig.add_subplot(2, 1, i+1) for i in range(2)]
axs = []
for _ax in _axs:
  divider = make_axes_locatable(_ax)
  axs.append(_ax)
  axs.append(divider.append_axes("bottom", size="20%", pad=0.0, sharex=_ax))
axs.append(divider.append_axes("left", size="20%", pad=0.0, sharey=axs[2]))
I = 5

for pft in ['CRO', 'DBF', 'GRA', 'ENF', 'CSH']:
  _df = df.loc[df.pft == pft, :]
  pft_leaf(_df, axs)

# df.groupby('pft').apply(pft_leaf, axs)
axs[1].set_xlabel('VPD (Pa)')
axs[3].set_xlabel('LAI')
axs[0].set_ylabel(paren_string)
plt.setp(axs[2].get_yticklabels(), visible=False)
axs[-1].set_ylabel(r'VPD')#$_{ETmin}$')
axs[2].set_title('VPD where ET = Min(ET) '\
                 r'($\frac{\partial \; ET}{\partial \; D} = 0$)')
axs[0].plot(axs[2].get_xlim(), [0., 0.], 'k--', linewidth=0.2)
# axs[2].set_ylim([0., np.around(df.vpd.quantile(q=0.95), decimals=-2)])
axs[2].set_ylim([0., 4000.])
axs[1].set_ylim([0.5,5.5])
axs[3].set_ylim([0.5,5.5])
axs[1].get_yaxis().set_visible(False)
axs[3].get_yaxis().set_visible(False)
axs[-1].set_xlim([0.5,5.5])
axs[-1].get_xaxis().set_visible(False)

axs[2].text(0.35, 250., r'$\frac{\partial \; ET}{\partial \; D} < 0$',\
            fontsize=15)
axs[2].text(1.4, 3500., r'$\frac{\partial \; ET}{\partial \; D} > 0$',\
            fontsize=15)
# axs[3].set_xlim(axs[2].get_xlim())
# axs[1].set_xlim(axs[0].get_xlim())

# axs[1].set_ylabel(r'-$\frac{\gamma c_s }{LAI \; 1.6 \; R\;  uWUE }$')
# axs[2].set_xlabel('VPD (Pa)')
# axs[3].set_xlabel('VPD (Pa)')
# axs[2].set_ylabel(paren_string)
# axs[3].set_ylabel(r'-$\frac{2 g_1 + \sqrt{D}}{2 (g_1 + \sqrt{D})^2}$')
axs[0].set_ylim((-1.8058955891452384, 0.87408219563044132))
for ax in axs[:1]:
  h, l = ax.get_legend_handles_labels()
  ax.legend(h, l, loc='best', fontsize=9)
# plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../doc/paper/fig05.pdf')


###### table 5 ####
importlib.reload(calc)
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
print('\n mean d_et *std vpd\n', mean_df['d_et_bar_std'])
print('\n mean d_et / rad\n', mean_df['d_et_bar_norm_rn'])
print('\n portion d_et < 0 \n',\
      counts)
mean_df['counts'] = counts
test = mean_df.loc[:, ['d_et', 'd_et_bar', 'd_et_bar_std',\
                       'd_et_bar_norm_rn', 'counts']]

print(test)

### Figure 6 ###
# note below really takes a long time t
def d_et_lai_fixed(_df):
  """returns d_et calced with mean lai"""
  temp_df = _df.copy()
  temp_df['lai'] = _df.lai.mean()
  # print(temp_df['lai'])
  _df['d_et_lai_fixed'] = calc.d_et(temp_df)
  temp_df['c_a'] = _df.c_a.mean()
  _df['d_et_lai_c_a_fixed'] = calc.d_et(temp_df)
  temp_df['gamma'] = _df.gamma.mean()
  _df['d_et_lai_all_fixed'] = calc.d_et(temp_df)
  # temp_df['gamma'] = _df.gamma.mean()
  # _df['d_et_lai_all_fixed'] = calc.d_et(temp_df)
  temp_df['c_a'] = _df.c_a
  _df['d_et_lai_gamma_fixed'] = calc.d_et(temp_df)
  return _df

df = df.groupby('pft').apply(d_et_lai_fixed)

def rh_d_et_min(_df):
  """returns rh where ET is minimized, as a function of T and LAI"""
  mean_df = _df.loc[:, ['r_moist', 'gamma', 'c_a',\
                        'lai', 'uwue_norm', 'g1']].mean()
  vpd = et_min_vpd(mean_df, mean_df.lai)
  t = np.linspace(_df.t_a.min(), _df.t_a.max())
  esat = met.vapor_pres(t)*100.
  print('mean esat', esat.mean())
  print('mean vpd', vpd)
  rh = (1. - vpd/esat)*100.
  rh[rh < 0.] = np.nan
  return t, rh

def make_ax_plot(_ax, var, _df, meta):
  """makes an axis plot"""
  # divider = make_axes_locatable(_ax)
  # axs.append(_ax)
  # _ax2 = divider.append_axes("right", size="20%", pad=0.0)

  vmax = meta['vmax']
  vmin = -vmax
  color = _ax.scatter(_df[meta['x_axis']], _df['t_a'], c=var, alpha=0.1,\
                      s=meta['size'], cmap=meta['cmap'],\
                      vmin=vmin, vmax=vmax)
  t, rh = rh_d_et_min(_df)
  # _ax.plot(rh, t, 'k-')
  if (meta['x_axis'] == 'vpd'):
    t_a = np.linspace(_df['t_a'].min(),_df['t_a'].max(), 200.)
    test = met.vapor_pres(t_a)*100.*(1. - 0.90)
    _ax.plot(test, t_a, 'k-')
    test = met.vapor_pres(t_a)*100.*(1. - 0.2)
    _ax.plot(test, t_a, 'k-')
  _ax.set_xlabel(meta['x_axis'])
  _ax.set_ylabel('T (C)')
  _ax.set_title('PFT: %s; %s'\
                % (str(_df['pft'][0]),\
                   meta['title']))
  cbar = plt.colorbar(color, ax=_ax)# , ax=_ax2)
  cbar.set_label(meta['title'])
  return

def scatter_plot_paper(_df, meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """

  #meta['x_axis'] = 'vpd'
  nplots = meta['nplots'] #5
  meta['size'] = 1
  meta['cmap'] = 'RdBu'

  fig = plt.figure()
  fig.set_figwidth(fig.get_figwidth()*nplots)

  if nplots == 4:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a']]# ,\
  elif nplots == 5:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}'\
              r'{g_a \partial \; D}(\overline{LAI, \gamma})$']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a'],\
             _df['d_et_lai_gamma_fixed']/_df['g_a']]
  else:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}'\
              r'{g_a \partial \; D}(\overline{LAI, \gamma})$',\
              'c_a_fixed',\
              'c_a and gamma fixed']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a'],\
             _df['d_et_lai_gamma_fixed']/_df['g_a'],\
             _df['d_et_lai_c_a_fixed']/_df['g_a'],\
             _df['d_et_lai_all_fixed']/_df['g_a']]

  axs = [fig.add_subplot(1, nplots, i+1) for i in range(nplots)]

  for ax, var, meta['title'] in zip(axs, _vars, titles):
    meta['vmax'] = np.nanmax(np.absolute([var.mean() +  2.*var.std(),\
                                          var.mean() - 2.*var.std()]))
    make_ax_plot(ax, var, _df, meta)

  plt.tight_layout()

  fname = '%s/climate_et/paper_plots/scatter/%s_%s.png'\
          % (os.environ['PLOTS'], _df.pft.iloc[0], meta['x_axis'])

  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir -p %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

plt.close('all')
meta = {}
meta['nplots'] = 4 # 5 4
meta['x_axis'] = 'vpd'
df.groupby('pft').apply(scatter_plot_paper, meta)
os.system('convert -append %s/climate_et/paper_plots/scatter/*%s.png '\
          '../doc/paper/fig06.png'\
          % (os.environ['PLOTS'], meta['x_axis']))


meta['x_axis'] = 'rh'
df.groupby('pft').apply(scatter_plot_paper, meta)
os.system('convert -append %s/climate_et/paper_plots/scatter/*_%s.png '\
          '../doc/paper/fig06b.png'\
          % (os.environ['PLOTS'], meta['x_axis']))


##### Table 5  #####
columns = []
for i in range(1, 5):
  columns.append('term_%d' % i)
  mean_df[columns[-1]] = 1./(mean_df.lai*mean_df.uwue_norm*mean_df.g1**i)

print(mean_df.loc[:, columns])
