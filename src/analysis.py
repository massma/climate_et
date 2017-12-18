import codebase.data_io as d_io
import os
import importlib
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

prep_df = d_prep.generate_vars(df)
df.to_pickle('%s/changjie/prepped_data.pkl' % os.environ['DATA'])
importlib.reload(d_calc)
dfs = d_calc.all_diagnostics(prep_df)
df.to_pickle('%s/changjie/diagnosed_data.pkl' % os.environ['DATA'])
importlib.reload(test)
test.run_all_tests(dfs)


