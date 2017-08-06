#! ~/edward/bin/python
"""
This script tests that changjie's method fits ours
"""
import time
import os
import glob
import importlib
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io as io

importlib.reload(pm)
WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
           comment='#', delimiter=',')
WUE.index = WUE.PFT
LV = 2.5e6

sitelist = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
sitelist.index = sitelist.Site

filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])
_ds = np.ones(len(filenames))*np.nan
_coef = pd.DataFrame(data={'PFT' : _ds, 'g0' : _ds, 'g1' : _ds,\
                           'count' : _ds, 'r2': _ds}, index=filenames)

time_start = time.time()
for i,filename in enumerate(filenames[:]):
  print('working on %d' % i)
  data = io.loadmat(filename)
  pft = str(np.squeeze(data['cover_type']))
  try:
    _ = WUE.loc[pft, :]
  except KeyError:
    print("file %s 's pft has no uWUE, moving on" % filename)
    continue

  atmos = {}
  atmos['t_a'] = np.squeeze(data['TA']) #C
  atmos['rh'] = np.squeeze(data['RH']*100.) #percent
  atmos['h'] = np.squeeze(data['H'])
  atmos['r_n'] = np.squeeze(data['NETRAD'])
  atmos['ustar'] = np.squeeze(data['USTAR']) #m/s
  atmos['p_a'] = np.squeeze(data['PA']*1000.) #Pa
  atmos['u_z'] = np.squeeze(data['WS'])

  canopy = {}
  canopy['g_flux'] = np.squeeze(data['G'])
  canopy['height'] = float(sitelist.loc[data['sitecode'], 'Canopy_h'])
  canopy['zmeas'] = float(sitelist.loc[data['sitecode'], 'Measure_h'])
  canopy['r_s'] = np.squeeze(data['Rs']) #s/m

  et = pm.penman_monteith(atmos, canopy)


  if np.isnan(np.nanmax(et)):
    print('error no et for %s' % filename)
  else:
    plt.figure()
    sns.regplot(x=data['LE'], y=et)
    plt.savefig('%s/%05d_garb.png' % (os.environ['PLOTS'], i))
    plt.show(block=False)

  atmos, canopy = pm.penman_monteith_prep(atmos, canopy)
  canopy['r_s'] = ((((atmos['delta']*(atmos['r_n']-canopy['g_flux'])+\
                      (1012.*atmos['rho_a']*atmos['vpd'])/atmos['r_a'])\
                     /np.squeeze(data['LE'])-atmos['delta'])\
                    /atmos['gamma'])-1.)*atmos['r_a']
  #print(data['Rs'][~np.isnan(data['Rs'])].shape)
  #print(canopy['r_s'][~np.isnan(canopy['r_s'])].shape)
  et = pm.penman_monteith(atmos, canopy)
  if np.isnan(np.nanmax(et)):
    print('error no et for %s' % filename)
  else:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sns.regplot(x=data['LE'], y=et, fit_reg=True, ax=ax)
    ax.plot(ax.get_xlim(), ax.get_xlim())
    plt.savefig('%s/%05d_garb_myfit.png' % (os.environ['PLOTS'], i))
    plt.show(block=False)
