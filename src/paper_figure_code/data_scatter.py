#! ~/edward/bin/python
"""
This script makes fig 6
"""
from shared_functions import *
import time

# mean_df is defined up at this level
local_fontsize = fontsize+3
### Figure 6 ###
# note below really takes a long time
def d_et_uwue_fixed(_df):
  """returns d_et calced with mean lai"""
  mean_row = mean_df.loc[_df.pft.iloc[0], :]
  temp_df = _df.copy()
  temp_df['uwue'] = mean_row.uwue
  _df['d_et_uwue_fixed'] = d_calc.sign(temp_df)*d_calc.scaling(temp_df)
  return _df

df = df.groupby('pft').apply(d_et_uwue_fixed)

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
  vmax = meta['vmax']
  vmin = -vmax
  _df['var'] = var
  _df = _df.assign(x_cut=pd.cut(_df['vpd'], meta['nbins']),\
                   t_a_cut=pd.cut(_df['t_a'], meta['nbins']))
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
  # print(_xstr)
  _ystr = grouped.index.levels[1].values
  # print(_ystr)
  # add some x, y modifies here to just grab the appropriate edge
  _x = np.array([interval.left for interval in _xstr])
  _xmax = np.array([interval.right for interval in _xstr]).max()
  _y = np.array([interval.left for interval in _ystr])
  _ymax = np.array([interval.right for interval in _ystr]).max()
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
  z = np.ma.masked_array(grouped.values, mask=np.isnan(grouped.values))
  color = _ax.pcolormesh(_x, _y, z, cmap=meta['cmap'],\
                         vmin=vmin, vmax=vmax)
  if str(_df['pft'].iloc[0]) == 'GRA':
    _ax.set_xlabel('%s (Pa)' % 'vpd'.upper(), fontsize=local_fontsize)
  if meta['title'] ==  r'$\frac{\partial \; ET}{\partial \, VPD}$':
    _ax.set_ylabel('T (C)', fontsize=local_fontsize)
  # _ax.set_title('%s: %s'\
  #               % (str(_df['pft'].iloc[0]),\
  #                  meta['title']), fontsize=local_fontsize+3)
  _ax.set_title('%s'\
                % (name_dict[_df['pft'].iloc[0]]), fontsize=local_fontsize+3)
  _ax.yaxis.set_tick_params(labelsize=local_fontsize-3)
  _ax.xaxis.set_tick_params(labelsize=local_fontsize-3)
  cbar = plt.colorbar(color, ax=_ax)# , ax=_ax2)
  cbar.set_label(meta['title'], fontsize=local_fontsize+3)
  y_lim = _ax.get_ylim()
  vpd_crit = et_min_vpd(mean_df.loc[_df.pft.iloc[0], :])
  _ax.plot(np.ones(2)*vpd_crit, np.array(y_lim), 'k-', linewidth=2.0)
  xlim = [0.0, dfs['95'].loc[_df.pft.iloc[0], 'vpd']+100.0]
  ylim = [dfs['5'].loc[_df.pft.iloc[0], 't_a']-1.0,\
          dfs['95'].loc[_df.pft.iloc[0], 't_a']+1.0]
  _ax.set_xlim(xlim)
  _ax.set_ylim(ylim)
  return

def scatter_plot_paper(_df, meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  nplots = 4
  meta['size'] = 1
  meta['cmap'] = 'RdBu'

  fig = plt.figure()
  fig.set_figwidth(fig.get_figwidth()*nplots)

  titles = [r'$\frac{\partial \; ET}{\partial \, VPD}$',\
            r'$\frac{\partial \; ET}{g_a \partial \, VPD}$',\
            r'$\frac{\partial \; ET}{\partial \, VPD}(\overline{\sigma})$',\
            r'$\frac{\partial \; ET}{g_a \partial \, VPD}'\
            r'(\overline{\sigma})$']
  _vars = [_df['d_et'],\
           _df['d_et']/_df['g_a'],\
           _df['d_et_uwue_fixed'],\
           _df['d_et_uwue_fixed']/_df['g_a']]# ,\

  axs = [fig.add_subplot(1, nplots, i+1) for i in range(nplots)]

  for i, ax, var, meta['title'] in zip(range(len(_vars)), axs, _vars, titles):
    ### set axis limits
    # if ((i % 2) == 0):
    #   j = 2
    # else:
    #   j = 3
    # meta['vmax'] = np.nanmax(np.absolute([_vars[j].mean()+2.*_vars[j].std(),\
    #                                       _vars[j].mean()-2.*_vars[j].std()]))
    meta['vmax'] = np.nanmax(np.absolute([var.mean()+1.5*var.std(),\
                                          var.mean()-1.5*var.std()]))
    make_ax_plot(ax, var, _df, meta)

  plt.tight_layout()

  fname = '%s/climate_et/paper_plots/scatter/%s_vpd.png'\
          % (os.environ['PLOTS'], _df.pft.iloc[0])

  try:
    plt.savefig(fname)
  except FileNotFoundError:
    os.system('mkdir -p %s' % '/'.join(fname.split('/')[:-1]))
    plt.savefig(fname)
  plt.show(block=False)
  return

start = time.time()

plt.close('all')
meta = {}
meta['sample'] = '' # 'sampled'


meta['nbins'] = 1000
# much faster if bins smaller, useful for testing
# meta['nbins'] = 10
df.groupby('pft').apply(scatter_plot_paper, meta)
os.system('convert -append %s/climate_et/paper_plots/scatter/*vpd.png '\
          '../../doc/paper/data_scatter%s.png'\
          % (os.environ['PLOTS'], meta['sample']))
os.system('cp ../../doc/paper/data_scatter%s.png '\
          '../../doc/paper/data_scatter%s.bak'
          % (meta['sample'], meta['sample']))

print('done with vpd, time was %f min' % ((time.time()-start)/60.))

