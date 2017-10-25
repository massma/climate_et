#! ~/edward/bin/python
"""
This script makes fig 6
"""
from shared_functions import *

### Figure 6 ###
# note below really takes a long time t
def d_et_lai_fixed(_df):
  """returns d_et calced with mean lai"""
  temp_df = _df.copy()
  temp_df['lai'] = _df.lai.mean()
  # print(temp_df['lai'])
  _df['d_et_lai_fixed'] = calc.d_et(temp_df)
  temp_df['c_a'] = _df.c_a.mean()
  _df['d_et_lai_c_a_fixed'] = calc.d_et(temp_df)
  temp_df['gamma'] = _df.gamma.mean()
  _df['d_et_lai_all_fixed'] = calc.d_et(temp_df)
  # temp_df['gamma'] = _df.gamma.mean()
  # _df['d_et_lai_all_fixed'] = calc.d_et(temp_df)
  temp_df['c_a'] = _df.c_a
  _df['d_et_lai_gamma_fixed'] = calc.d_et(temp_df)
  return _df

df = df.groupby('pft').apply(d_et_lai_fixed)

def rh_d_et_min(_df):
  """returns rh where ET is minimized, as a function of T and LAI"""
  mean_df = _df.loc[:, ['r_moist', 'gamma', 'c_a',\
                        'lai', 'uwue_norm', 'g1']].mean()
  vpd = et_min_vpd(mean_df, mean_df.lai)
  t = np.linspace(_df.t_a.min(), _df.t_a.max())
  esat = met.vapor_pres(t)*100.
  print('mean esat', esat.mean())
  print('mean vpd', vpd)
  rh = (1. - vpd/esat)*100.
  rh[rh < 0.] = np.nan
  return t, rh

def meshgrid_apply(_df, column='var', sample=False):
  """calculates a mean on column if size > 1"""
  #if _df.size > 1:
  if sample:
    # print('sampled!')
    return float(_df.loc[:, column].sample(n=1))
  else:
    return _df.loc[:, column].mean(axis=0)

def make_ax_plot(_ax, var, _df, meta):
  """makes an axis plot"""
  # divider = make_axes_locatable(_ax)
  # axs.append(_ax)
  # _ax2 = divider.append_axes("right", size="20%", pad=0.0)

  vmax = meta['vmax']
  vmin = -vmax
  # color = _ax.scatter(_df[meta['x_axis']], _df['t_a'], c=var, alpha=0.1,\
  #                     s=meta['size'], cmap=meta['cmap'],\
  #                     vmin=vmin, vmax=vmax)
  #below is a hack, should have done by name to begin with
  _df['var'] = var
  _df = _df.assign(x_cut=pd.cut(_df[meta['x_axis']], 1000),\
                   t_a_cut=pd.cut(_df['t_a'], 1000))
  if meta['sample'] == 'sampled':
    grouped = _df.groupby(['x_cut', 't_a_cut']).apply(meshgrid_apply,\
                                                      sample=True)
  else:
    grouped = _df.groupby(['x_cut', 't_a_cut']).apply(meshgrid_apply)
  grouped = grouped.reset_index()
  grouped.columns = ['x_cut', 't_a_cut', 'var']
  grouped = grouped.pivot('x_cut', 't_a_cut').transpose()
  # print(grouped)
  _xstr = grouped.columns.values
  _ystr = grouped.index.levels[1].values
  # add some x, y modifies here to just grab the appropriate edge
  _x = np.array([float(s[1:-1].split(', ')[0]) for s in _xstr])
  _xmax = np.array([float(s[1:-1].split(', ')[-1]) for s in _xstr]).max()
  _y = np.array([float(s[1:-1].split(', ')[0]) for s in _ystr])
  _ymax = np.array([float(s[1:-1].split(', ')[-1]) for s in _ystr]).max()
  grouped.columns = _x
  grouped.index = _y
  grouped = grouped.sort_index(axis=1)
  # print('grouped sortex axis 1', grouped)
  grouped = grouped.sort_index(axis=0)
  # print('grouped sort axis 2', grouped)
  _x = grouped.columns.values
  _y = grouped.index.values
  _x = np.append(_x, _xmax)
  _y = np.append(_y, _ymax)
  _x, _y = np.meshgrid(_x, _y)
  # grouped[np.isnan(grouped)] = 0.
  # print('_x mesh', _x)
  # print('_y mesh', _y)
  z = np.ma.masked_array(grouped.values, mask=np.isnan(grouped.values))
  color = _ax.pcolormesh(_x, _y, z, cmap=meta['cmap'],\
                         vmin=vmin, vmax=vmax)
  t, rh = rh_d_et_min(_df)
  # _ax.plot(rh, t, 'k-')n
  if (meta['x_axis'] == 'vpd'):
    t_a = np.linspace(_df['t_a'].min(),_df['t_a'].max(), 200.)
    test = met.vapor_pres(t_a)*100.*(1. - 0.90)
    _ax.plot(test, t_a, 'k-')
    test = met.vapor_pres(t_a)*100.*(1. - 0.2)
    _ax.plot(test, t_a, 'k-')
  _ax.set_xlabel(meta['x_axis'])
  _ax.set_ylabel('T (C)')
  _ax.set_title('PFT: %s; %s'\
                % (str(_df['pft'][0]),\
                   meta['title']))
  cbar = plt.colorbar(color, ax=_ax)# , ax=_ax2)
  cbar.set_label(meta['title'])
  return

def scatter_plot_paper(_df, meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """

  #meta['x_axis'] = 'vpd'
  nplots = meta['nplots'] #5
  meta['size'] = 1
  meta['cmap'] = 'RdBu'

  fig = plt.figure()
  fig.set_figwidth(fig.get_figwidth()*nplots)

  if nplots == 4:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a']]# ,\
  elif nplots == 5:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}'\
              r'{g_a \partial \; D}(\overline{LAI, \gamma})$']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a'],\
             _df['d_et_lai_gamma_fixed']/_df['g_a']]
  else:
    titles = [r'$\frac{\partial \; ET}{\partial \; D}$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}$',\
              r'$\frac{\partial \; ET}{\partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \; ET}{g_a \partial \; D}(\overline{LAI})$',\
              r'$\frac{\partial \;n ET}'\
              r'{g_a \partial \; D}(\overline{LAI, \gamma})$',\
              'c_a_fixed',\
              'c_a and gamma fixed']
    _vars = [_df['d_et'],\
             _df['d_et']/_df['g_a'],\
             _df['d_et_lai_fixed'],\
             _df['d_et_lai_fixed']/_df['g_a'],\
             _df['d_et_lai_gamma_fixed']/_df['g_a'],\
             _df['d_et_lai_c_a_fixed']/_df['g_a'],\
             _df['d_et_lai_all_fixed']/_df['g_a']]

  axs = [fig.add_subplot(1, nplots, i+1) for i in range(nplots)]

  for ax, var, meta['title'] in zip(axs, _vars, titles):
    meta['vmax'] = np.nanmax(np.absolute([var.mean() +  2.*var.std(),\
                                          var.mean() - 2.*var.std()]))
    make_ax_plot(ax, var, _df, meta)

  plt.tight_layout()

  fname = '%s/climate_et/paper_plots/scatter/%s_%s.png'\
          % (os.environ['PLOTS'], _df.pft.iloc[0], meta['x_axis'])

  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir -p %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

plt.close('all')
meta = {}
meta['nplots'] = 4 # 5 4
meta['x_axis'] = 'rh'
meta['sample'] = 'sampled'
# meta['sample'] = ''
df.groupby('pft').apply(scatter_plot_paper, meta)
os.system('convert -append %s/climate_et/paper_plots/scatter/*_%s.png '\
          '../../doc/paper/fig06b%s.png'\
          % (os.environ['PLOTS'], meta['x_axis'], meta['sample']))


meta['x_axis'] = 'vpd'
df.groupby('pft').apply(scatter_plot_paper, meta)
os.system('convert -append %s/climate_et/paper_plots/scatter/*%s.png '\
          '../../doc/paper/fig06%s.png'\
          % (os.environ['PLOTS'], meta['x_axis'], meta['sample']))

