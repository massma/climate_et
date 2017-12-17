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
importlib.reload(d_calc)
dfs = d_calc.all_diagnostics(prep_df)
importlib.reload(test)
test.run_all_tests(dfs)


# importlib.reload(test)
# test.test_et_model(calc_df)



