import codebase.data_io as d_io
import importlib

importlib.reload(d_io)

df = d_io.load_mat_data()
# prepped_data = prep_data(data)
# analyzed_data = analyze_data(prepped_data)

