import codebase.data_io as d_io
import os
import importlib
import pickle
import numpy as np
import pandas as pd
import codebase.data_prep as d_prep
import codebase.data_calc as d_calc
import codebase.test_functions as test

importlib.reload(d_io)
importlib.reload(d_prep)
importlib.reload(d_calc)

# df = d_io.load_mat_data()
# df.to_pickle('%s/changjie/loaded_data.pkl' % os.environ['DATA'])
df = pd.read_pickle('%s/changjie/loaded_data.pkl' % os.environ['DATA'])
print('pre-prep data shape:', df.shape)

prep_df = d_prep.generate_vars(df)
print('pre-calc data shape:', prep_df.shape)
df.to_pickle('%s/changjie/prepped_data.pkl' % os.environ['DATA'])
importlib.reload(d_calc)
dfs = d_calc.all_diagnostics(prep_df)
prep_df = d_prep.generate_vars(df)
print('post-calc data shape:', dfs['full'].shape)

with open('%s/changjie/diagnosed_data.pkl' % os.environ['DATA'],\
          mode='wb') as file :
  pickle.dump(dfs, file)
importlib.reload(test)
test.run_all_tests(dfs)
