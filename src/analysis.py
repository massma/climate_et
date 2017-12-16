import codebase.data_io as d_io
import importlib
import codebase.data_prep as d_prep

importlib.reload(d_io)


# df = d_io.load_mat_data()
# df.to_pickle('%s/changjie/loaded_data.pkl')

importlib.reload(d_prep)
prep_df = d_prep.generate_vars(df)
# analyzed_data = analyze_data(prepped_data)

