#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import time
import os
import importlib
import pandas as pd
import codebase.penman_monteith as pm
import codebase.data_io as d_io

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

def main():
  """wrapper for main script"""
  start = time.time()
  coef = pd.read_csv('../dat/site_coef_mm_s_medlyn.csv')
  coef.index = coef['Unnamed: 0']
  outdir = '%s/changjie/pandas_data/' % os.environ['DATA']
  for i, index in enumerate(coef.index[:]):
    print('Working on %s, file number %d, time elapsed: %f m' \
          % (index, i, (time.time()-start)/60.))
    atmos, canopy, data = d_io.load_mat_data(index)
    canopy['stomatal_model'] = 'adam_medlyn'
    # belof converts from mm/s to m/s, only do g0 b/c of funcitonal form
    canopy['g0'] = coef.loc[index, 'g0']/1000.
    canopy['g1'] = coef.loc[index, 'g1']
    data['et'] = pm.recursive_penman_monteith(atmos, canopy)
    atmos['vpd_leaf'] = atmos['vpd'] + 1.0
    data['et_leaf'] = pm.recursive_penman_monteith(atmos, canopy)
    atmos['vpd_leaf'] = atmos['vpd']
    atmos['vpd'] = atmos['vpd'] + 1.0
    data['et_atm'] = pm.recursive_penman_monteith(atmos, canopy)
    atmos['vpd_leaf'] = atmos['vpd_leaf'] + 1.0
    data['et_all'] = pm.recursive_penman_monteith(atmos, canopy)
    dfout = pd.concat([atmos, canopy, data], axis=1)
    fname = ''.join(index.split('/')[-1].split('.')[:-1])
    dfout.to_pickle('%s/%s.pkl' % (outdir, fname))
  print('time was %f s' % ((time.time()-start)))
  return

if str(__name__) == '__main__':
  main()
