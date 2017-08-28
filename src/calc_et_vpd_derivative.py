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

def leaf_vpd(atmos, canopy, lai):
  """calculates the leaf term in dET/dVPD (see doc folder)"""
  return atmos['gamma']*atmos['c_a']*atmos['p_a']*/\
    (lai*1.6*pm.R_STAR*(273.15+atmos['t_a'])*canopy['uwue'])\
    *(2.*canopy['g1'] + np.sqrt(atmos['vpd']))\
    /(2.*(canopy['g1'] + np.sqrt(atmos['vpd']))**2)


def calc_derivative(atmos, canopy, data):
  """adds various derivative fields to data, given atmos and canopy"""
  data['et'] = pm.penman_monteith_uwue(atmos, canopy)
  data['scaling'] = 1./(atmos['r_a']*(atmos['delta'] + atmos['gamma']))
  data['vpd_atm'] = atmos['rho_a']*pm.CP
  data['vpd_leaf'] = leaf_vpd(atmos, canopy, canopy['lai'])
  data['lai_fit'] = atmos['gamma']*atmos['c_a']\
                    *np.sqrt(atmos['vpd'])*atmos['p_a']\
                    /(atmos['r_a']\
                      *(data['et_obs']*(atmos['gamma']+atmos['delta'])\
                        -atmos['delta']*(atmos['r_n']-canopy['g_flux'])\
                        -1./atmos['r_a']*atmos['rho_a']*pm.CP*atmos['vpd'])
                      *1.6*pm.R_STAR*(273.15 + canopy['t_a'])\
                      *canopy['uwue']*(1.+canopy['g1']/atmos['vpd']))
  data['vpd_leaf_hourly'] = leaf_vpd(atmos, canopy, data['lai_fit'])
  return data

#def main():
"""wrapper for main script"""


start = time.time()
coef = pd.read_csv('../dat/site_coef_mm_s_W_m2_medlyn_lai_nolim.csv')
coef.index = coef['Unnamed: 0']
outdir = '%s/changjie/pandas_data_lai_fit_nolim/' % os.environ['DATA']
for i, index in enumerate(coef.index[:1]):
  print('Working on %s, file number %d, time elapsed: %f m' \
        % (index, i, (time.time()-start)/60.))
  atmos, canopy, data = d_io.load_mat_data(index)
  canopy['et_stomatal_model'] = 'adam_medlyn'
  canopy['lai'] = coef.loc[index, 'lai']
  canopy['g1'] = coef.loc[index, 'g1']
  canopy['uwue'] = pm.WUE.loc[canopy['pft'].iloc[0], 'u_wue_yearly']
  data = calc_derivative(atmos, canopy, data)
  dfout = pd.concat([atmos, canopy, data], axis=1)
  fname = ''.join(index.split('/')[-1].split('.')[:-1])
  dfout.to_pickle('%s/%s.pkl' % (outdir, fname))
print('time was %f s' % ((time.time()-start)))
#return

# if str(__name__) == '__main__':
#   main()
