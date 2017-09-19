#! ~/edward/bin/python
"""
This script makes some simple analytical plots of vars that
vary in the dET/D_s term
"""
import glob
import os
import pandas as pd
import numpy as np
import metcalcs as met
import seaborn as sns
import resource
from scipy.stats import spearmanr
nfrom scipy.optimize import least_squares
import matplotlib.pyplot as plt
import matplotlib as mpl
import codebase.penman_monteith as pm
import util
import tensorflow as tf
import edward as ed
from edward.models import Normal
from sklearn.model_selection import train_test_split
from datetime import datetime
import time
mpl.rcParams.update(mpl.rcParamsDefault)

df = pd.read_pickle('%s/changjie/full_pandas_lai_clean.pkl'\
                    % os.environ['DATA'])

df['g_a'] = 1./df['r_a']

_time = pd.DatetimeIndex(df.time)
df['hour'] = _time.hour
df['jd'] = _time.dayofyear
df['year'] = _time.year

def get_uwue(_df):
  """lazy way to get uwue"""
  return _df.loc[:, ['uwue']].iloc[0]

def d_et(_df):
  """calcs the full d ET/d Ds, confirmed correct vs df"""
  return _df['g_a']*_df['p_a']/\
    ((_df['t_a']+ 273.15)*(_df['gamma']+_df['delta']))*\
    (pm.CP/_df['r_moist']-_df['gamma']*_df['c_a']*pm.LV/\
     (_df['lai']*1.6*pm.R_STAR*_df['uwue'])*\
     (2.*_df['g1']+np.sqrt(_df['vpd']))\
     /(2.*(_df['g1']+np.sqrt(_df['vpd']))**2))

FIG = plt.figure()

def fit_curves(_df, ax):
  """fits curves"""
  jd = _df.jd.values
  lai = _df.lai.values
  def seasonal_lai(x):
    """
    model to fit to lai,
    x: x[0] * sin(jd/365*2.pi + x[1]) + x[2]
    """
    return x[0]*np.sin(jd/365.*2.*np.pi + x[1]) + x[2] - lai
  x0 = [lai.std(), 0., lai.mean()]
  bounds = ([0., 0., 0.], [10., 2.*np.pi, 10.])
  result = least_squares(seasonal_lai, x0, bounds=bounds)
  if result['success']:
    seasonal = seasonal_lai(result['x']) + lai
    df_out = pd.DataFrame(data={'modelled_cycle': seasonal}, index=_df.index)
    for i, _x in enumerate(result['x']):
      df_out['x%d' % i] = _x
    ax.plot(_df.time.values, seasonal, 'k-')
    return df_out
  else:
    print('error, no convergent solution for %s at year %d'\
          % (_df.site.iloc[0], int(_df.year.iloc[0])))
    return np.ones(3)*np.nan
  return

def plot_lai(_df, fit_plot=False, suff=''):
  ax = FIG.add_subplot(111)
  ax.scatter(_df.time.values, _df.lai.values, s=1)
  FIG.autofmt_xdate()
  x = _df.groupby('year').apply(fit_curves, ax)
  util.test_savefig('%s/climate_et/lai_plots/%s_%s%s.png'\
                    % (os.environ['PLOTS'], _df.pft.iloc[0],\
                       _df.site.iloc[0], suff))
  FIG.clf()
  return x

# # below is to generate datta
# x = df.groupby('site').apply(plot_lai, fit_plot=True)
# full_df = pd.concat([df, x], axis=1)
# full_df.to_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
#                   % os.environ['DATA'])

df = pd.read_pickle('%s/changjie/full_pandas_seasonal_fit.pkl'\
                    % os.environ['DATA'])

# bestfit = df.loc[df.x0.argmax()]

# df['residual'] = df.lai-df.modelled_cycle

# def residual_plot(_df):
#   """plots residuals"""
#   ax = FIG.add_subplot(111)
#   _df['residual'].hist(ax=ax)
#   util.test_savefig('%s/climate_et/lai_plots/histogram/%s_%s_%d.png'\
#                     % (os.environ['PLOTS'], _df.pft.iloc[0],\
#                        _df.site.iloc[0], int(_df.year.iloc[0])))
#   FIG.clf()
#   return

# def year_group(_df):
#   """takes groupedby site and groups by year"""
#   _df.groupby('year').apply(residual_plot)
#   return

# df.groupby('site').apply(year_group)


columns = ['lai', 'jd', 'year',\
           'modelled_cycle', 'x0', 'x1', 'x2', 'site', 'pft']
subset = df.loc[(df.year == 2003) & (df.site == 'US-MMS'), columns]

subset['residual'] = subset.lai-subset.modelled_cycle
# well-fit subset to test with edward approahc
# plot_lai(subset, suff='_subset')

def gen_data(_df):
  """get testing and training data from _df"""
  train, test = train_test_split(_df, test_size=0.2)
  y_train = train.lai.values
  y_test = test.lai.values
  x_train = train.jd.values
  x_test = test.jd.values
  return x_train, x_test, y_train, y_test

x_train, x_test, y_train, y_test = gen_data(subset)

N = y_train.size

#initialize graph

# remove default graph
tf.reset_default_graph()

#model params
y_scale = subset.lai.std()
amp = Normal(loc=tf.ones(1)*subset.lai.std(), scale=tf.ones(1))
phase = Normal(loc=tf.ones(1)*np.pi, scale=tf.ones(1)*np.pi)
offset = Normal(loc=tf.ones(1)*subset.lai.mean(), scale=tf.ones(1))
noise = Normal(loc=tf.zeros(1), scale=tf.ones(1))
x = tf.placeholder(tf.float32, [N])
y = Normal(loc=amp*tf.sin(x/365.*2.*np.pi+phase) + offset,\
           scale=tf.nn.softplus(noise))

q_amp = Normal(loc=tf.Variable(tf.random_normal([1])),
                scale=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
q_phase = Normal(loc=tf.Variable(tf.random_normal([1])),
                 scale=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
q_offset = Normal(loc=tf.Variable(tf.random_normal([1])),
                  scale=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
q_noise = Normal(loc=tf.Variable(tf.random_normal([1])),
                  scale=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))

inference = ed.KLqp({amp: q_amp, phase: q_phase,\
                     offset : q_offset, noise: q_noise},\
                    data={x: x_train, y: y_train})
inference.run(n_samples=10, n_iter=1000)

y_post = ed.copy(y, {amp: q_amp, phase: q_amp, offset: q_offset})

print("Mean squared error on test data:")
print(ed.evaluate('mean_squared_error', data={x: x_train, y_post: y_train}))

print("Mean absolute error on test data:")
print(ed.evaluate('mean_absolute_error', data={x: x_train, y_post: y_train}))

def visualize(x_data, y_data, amp, phase, offset, n_samples=10, name=''):
  amp_samples = amp.sample(n_samples).eval()
  phase_samples = phase.sample(n_samples).eval()
  offset_samples = offset.sample(n_samples).eval()
  ax = FIG.add_subplot(111)
  ax.scatter(x_data, y_data)
  inputs = np.linspace(0., 365., num=400)
  ax.plot(inputs, subset.x0.iloc[0]*np.sin(inputs/365.*2.*np.pi\
                                           +subset.x1.iloc[0])\
          +subset.x2.iloc[0], 'k-', label='least-sq')
  ax.plot(subset.jd.values, subset.modelled_cycle.values,\
          'ko', label='least-sq')
  for ns in range(n_samples):
    output = amp_samples[ns]*np.sin(inputs/365.*2.*np.pi+phase_samples[ns])\
             + offset_samples[ns]
    ax.plot(inputs, output)
  plt.legend(loc='best')
  util.test_savefig('%s/climate_et/lai_plots/bayesian/%s_%d.png'\
                    % (os.environ['PLOTS'], name, int(time.time())))
  FIG.clf()
  return

visualize(x_train, y_train, amp, phase, offset, name='prior')
visualize(x_train, y_train, q_amp, q_phase, q_offset, name='posterior')

print('bayesian residual std: %f' % tf.nn.softplus(q_noise.loc).eval())
print('actual residual std: %f' % subset.residual.std())
