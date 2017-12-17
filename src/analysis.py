import codebase.data_io as d_io
import os
import importlib
import numpy as np
import codebase.data_prep as d_prep
import codebase.data_calc as d_calc
import codebase.test_functions as test

importlib.reload(d_io)
importlib.reload(d_prep)
importlib.reload(d_calc)

df = d_io.load_mat_data()
df.to_pickle('%s/changjie/loaded_data.pkl' % os.environ['DATA'])

prep_df = d_prep.generate_vars(df)
calc_df = d_calc.uwue(prep_df)

def max_diff(quant1, quant2):
  """caclualte the maximuma boslute difference"""
  return np.nanmax(np.absolute(quant1-quant2))

importlib.reload(calc_df)
test.test_et_model(calc_df)
#is uwue goo
