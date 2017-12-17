import codebase.data_io as d_io
import os
import importlib
import codebase.data_prep as d_prep
import codebase.data_calc as d_calc

importlib.reload(d_io)

df = d_io.load_mat_data()
df.to_pickle('%s/changjie/loaded_data.pkl' % os.environ['DATA'])

importlib.reload(d_prep)
prep_df = d_prep.generate_vars(df)
analyzed_data = data_calc.(prepped_data)

