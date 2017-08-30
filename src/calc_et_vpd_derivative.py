#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import glob
import time
import os
import importlib
import pandas as pd
import codebase.penman_monteith as pm
import codebase.data_io as d_io
import numpy as np

importlib.reload(d_io)
importlib.reload(pm)

WUE = pd.read_csv('../dat/zhou_et_al_table_4.csv',\
           comment='#', delimiter=',')
WUE.index = WUE.PFT
LV = 2.5e6

SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

def leaf_vpd(atmos, canopy, lai):
  """calculates the leaf term in dET/dVPD (see doc folder)"""
  return -atmos['gamma']*atmos['c_a']*pm.LV*atmos['p_a']/\
    (lai*1.6*pm.R_STAR*(273.15+atmos['t_a'])*canopy['uwue'])\
    *(2.*canopy['g1'] + np.sqrt(atmos['vpd']))\
    /(2.*(canopy['g1'] + np.sqrt(atmos['vpd']))**2)


def calc_derivative(atmos, canopy, data):
  """adds various derivative fields to data, given atmos and canopy"""
  canopy['lai'] = - atmos['gamma']*atmos['c_a']*pm.LV\
                *np.sqrt(atmos['vpd'])*atmos['p_a']\
                /(atmos['r_a']\
                  *(data['et_obs']*(atmos['gamma']+atmos['delta'])\
                    -atmos['delta']*(atmos['r_n']-canopy['g_flux'])\
                    -1./atmos['r_a']*atmos['rho_a']*pm.CP*atmos['vpd'])
                  *1.6*pm.R_STAR*(273.15 + atmos['t_a'])\
                  *canopy['uwue']*(1.+canopy['g1']/np.sqrt(atmos['vpd'])))
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  data['scaling'] = 1./(atmos['r_a']*(atmos['delta'] + atmos['gamma']))
  data['vpd_atm'] = atmos['rho_a']*pm.CP
  data['vpd_leaf'] = leaf_vpd(atmos, canopy, canopy['lai'])
  atmos['vpd'] = atmos['vpd'] + 1.0
  data['et_all'] = pm.penman_monteith_uwue(atmos, canopy)
  atmos['vpd'] = atmos['vpd'] - 1.0
  return data

#def main():
"""wrapper for main script"""


start = time.time()
outdir = '%s/changjie/pandas_data_lai/' % os.environ['DATA']

filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])

time_start = time.time()
for filename in filenames[:]:
  print('working on %s' % filename)
  atmos, canopy, data = d_io.load_mat_data(filename)
  if (data.et_obs.count() > 0) & (canopy.dropna().uwue.count() > 0):
    data = calc_derivative(atmos, canopy, data)
    dfout = pd.concat([atmos, canopy, data], axis=1)
    fname = ''.join(filename.split('/')[-1].split('.')[:-1])
    dfout.to_pickle('%s/%s.pkl' % (outdir, fname))
  else:
    print('filename %s is invalid, et count %d and canopy count %d'\
         % (filename, data.et_obs.count(), canopy.dropna().uwue.count()))

print('time was %f s' % ((time.time()-start)))
#return

# if str(__name__) == '__main__':
#   main()
