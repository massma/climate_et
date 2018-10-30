#! ~/edward/bin/python
"""
This script makes fig 6
"""
from shared_functions import *
from pdb import set_trace as bp
import time
mpl.rcParams.update(small_ax_params)
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
  grouped = grouped.reset_index(level=0, drop=True)
  _xstr = grouped.columns.values
  _ystr = grouped.index.values
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
  if _df.pft.iloc[0] == 'DBF':
    yticks = np.arange(10, 31, 5)
    _ax.set_yticks(yticks)
    _ax.set_yticklabels(['%2.0f' % y for y in yticks])
  custom_xlabel(_df, _ax, '%s (Pa)' % 'vpd'.upper(), fontsize=small_ax_fontsize)
  custom_ylabel(_df, _ax, 'T (C)', fontsize=small_ax_fontsize)
  _ax.set_title('%s'\
                % (name_dict[_df['pft'].iloc[0]]), fontsize=local_fontsize)
  # _ax.yaxis.set_tick_params(labelsize=local_fontsize-3)
  # _ax.xaxis.set_tick_params(labelsize=local_fontsize-3)
  cbar = plt.colorbar(color, ax=_ax)# , ax=_ax2)
  if ((_df.pft.iloc[0] == 'MF') | (_df.pft.iloc[0] == 'WSA')\
      | (_df.pft.iloc[0] == 'SAV')):
    cbar.set_label(meta['title'], fontsize=small_ax_fontsize)
  y_lim = _ax.get_ylim()
  # vpd_crit = et_min_vpd(mean_df.loc[_df.pft.iloc[0], :])
  # _ax.plot(np.ones(2)*vpd_crit, np.array(y_lim), 'k-', linewidth=2.0)
  xlim = [0.0, dfs['95'].loc[_df.pft.iloc[0], 'vpd']+100.0]
  ylim = [dfs['5'].loc[_df.pft.iloc[0], 't_a']-1.0,\
          dfs['95'].loc[_df.pft.iloc[0], 't_a']+1.0]
  _ax.set_xlim(xlim)
  _ax.set_ylim(ylim)
  return

def scatter_plot_paper(_df, ax, meta):
  """
  creates scatter of derivatives wrt to VPD, assumes Delta(vpd) = 1.0 Pa
  """
  meta['size'] = 1
  meta['cmap'] = 'RdBu'
  meta['title'] = r'$\frac{\partial \; ET}{\partial \, VPD}$ '\
                  'W m$^{-2}$ Pa$^{-1}$'

  var = _df['d_et']

  meta['vmax'] = np.nanmax(np.absolute([var.mean()+1.0*var.std(),\
                                          var.mean()-1.0*var.std()]))
  make_ax_plot(ax, var, _df, meta)

  return

start = time.time()

plt.close('all')
meta = {}
meta['sample'] = '' # 'sampled'

# much faster if bins smaller, useful for testing
meta['nbins'] = 1000

panel_wrapper(df, scatter_plot_paper, "data_scatter.png", args=(meta,))

print('done with vpd, time was %f min' % ((time.time()-start)/60.))
