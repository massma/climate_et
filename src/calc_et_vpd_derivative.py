#! ~/edward/bin/python
"""
This module uses site specific medlyn fits to calcualte
simulated ET, VPD and d ET/ d VPD for full, leaf and atm
"""
import time
import os
import glob
import importlib
import pandas as pd
import numpy as np
import codebase.penman_monteith as pm
from scipy.optimize import leastsq
import scipy.io as io
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

#def main():

"""wrapper for main script"""
coef = pd.read_csv('../dat/site_coef_mm_s_medlyn.csv')
coef.index = coef['Unnamed: 0']
for index in coef.index[2:3]:
  atmos, canopy, data = d_io.load_mat_data(index)
  canopy['stomatal_model'] = 'medlyn'
  # belof converts from mm/s to m/s, only do g0 b/c of funcitonal form
  canopy['g0'] = coef.loc[index, 'g0']/1000.
  canopy['g1'] = coef.loc[index, 'g1']
  canopy = canopy.iloc[:5,:]
  atmos = atmos.iloc[:5,:]
  et = pm.recursive_penman_monteith(atmos, canopy)
# if str(__name__) == '__main__':
#   coef = main()
