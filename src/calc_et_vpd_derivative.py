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
import codebase.calc_tools as calc
import numpy as np

importlib.reload(d_io)
importlib.reload(pm)

SITELIST = pd.read_csv('%s/changjie/fluxnet_algorithm/'\
                       'Site_list_(canopy_height).csv' % os.environ['DATA'],\
                       delimiter=',')
SITELIST.index = SITELIST.Site

H2O = 18.01528e-3 #molecular mass kg/mol


#def main():
"""wrapper for main script"""


start = time.time()
outdir = '%s/changjie/pandas_data_c3_c4/' % os.environ['DATA']

filenames = glob.glob('%s/changjie/MAT_DATA/*.mat' % os.environ['DATA'])

time_start = time.time()
for filename in filenames[:]:
  print('working on %s' % filename)
  atmos, canopy, data = d_io.load_mat_data(filename)
  if (data.et_obs.count() > 0) & (canopy.dropna().uwue.count() > 0):
    atmos, canopy, data = calc.calc_derivative(atmos, canopy, data)
    idx = ((canopy['lai'] < canopy['lai'].quantile(q=0.95)) & \
           (canopy['lai'] > canopy['lai'].quantile(q=0.05)))
    lai_std = canopy.loc[idx, 'lai'].std()/canopy.loc[idx, 'lai'].mean()
    if canopy['pft'].iloc[0] == 'GRA':
      canopy['g1'] = pm.WUE_MEDLYN.loc['C4G', 'g1M']
      _atmos, _canopy, _data = calc.calc_derivative(atmos, canopy, data)
      idx = ((_canopy['lai'] < _canopy['lai'].quantile(q=0.95)) & \
             (_canopy['lai'] > _canopy['lai'].quantile(q=0.05)))
      _lai_std = _canopy.loc[idx, 'lai'].std()/_canopy.loc[idx, 'lai'].mean()
      if _lai_std < lai_std:
        print('for filename %s c4 grass std (%f) is less than'\
              'c3 gras std (%f)' % (filename, _lai_std, lai_std))
        print('old d_et: %f' % data.d_et.mean())
        data = _data
        canopy = _canopy
        atmos = _atmos
        canopy.pft = 'C4G'
        print('new d_et: %f' % data.d_et.mean())
        print(canopy.g1.iloc[0])
      else:
        print('for filename %s c4 gras std (%f)'\
              'is greater than c3 grass std (%f)'\
              % (filename, _lai_std, lai_std))
        canopy['g1'] = pm.WUE_MEDLYN.loc['GRA', 'g1M']
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
# now concat and clean up LAI outliers
def concat_dfs(folder='pandas_data_v2', fname='full_pandas_v2'):
  """
  puts all the individual site data into one pdf, and adds a site column to df
  """
  dfs = []
  filenames = glob.glob('%s/changjie/%s/*'\
                        % (os.environ['DATA'], folder))
  for filename in filenames[:]:
    _df = pd.read_pickle(filename)
    _df['site'] = ''.join(filename.split('/')[-1].split('.')[:-1])
    dfs.append(_df)
  full_df = pd.concat(dfs)
  full_df.to_pickle('%s/changjie/%s.pkl'\
                    % (os.environ['DATA'], fname))
  return full_df

def clean_df(_df, var='lai'):
  """remove unphysical LAI values from a df"""
  out = _df.loc[((_df[var] > 0.1) & (_df[var] < 100.)), :]
  return out

def site_clean(_df, var='lai'):
  """this will remove some percentile of data"""
  out = _df.loc[((_df[var] < _df[var].quantile(q=0.95)) & \
                (_df[var] > _df[var].quantile(q=0.05))), :]
  return out

def print_cv_lai_lai_gpp(_df):
  """prints the cv of lai and lai_gpp"""
  print(_df.site.iloc[0], _df.pft.iloc[0])
  print('cv lai: %f' % (_df.lai.std()/_df.lai.mean()))
  print('cv lai_gpp: %f'  %(_df.lai_gpp.std()/_df.lai_gpp.mean()))
  return (_df.lai.std()/_df.lai.mean())# ,\
          # (_df.lai_gpp.std()/_df.lai_gpp.mean())
reload_data = True
if reload_data:
  concat_dfs(folder='pandas_data_c3_c4',\
             fname='full_pandas_c3_c4')
  df = pd.read_pickle('%s/changjie/full_pandas_c3_c4.pkl'\
                      % os.environ['DATA'])
  meta = {}
  meta['folder_label'] = 'site'
  meta['folder'] = 'hist_plots'
  meta['var'] = 'lai_gpp'
  print(df.shape)
  df = df.groupby('site').apply(site_clean)
  print(df.shape)
  df = clean_df(df)
  df = clean_df(df, var='lai_gpp')
  # test = df.groupby('site').apply(site_clean, 'lai_gpp')
  # test = clean_df(test, var='lai_gpp')
  df.to_pickle('%s/changjie/full_pandas_c3_c4_clean.pkl'\
               % os.environ['DATA'])

  cvs = df.groupby('site').apply(print_cv_lai_lai_gpp)
  print(df.shape)
  #df.groupby('site').apply(histogram, meta)
  # histogram(df, meta)
  # meta['var'] = 'lai'
  # histogram(df, meta)
