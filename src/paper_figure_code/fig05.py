#! ~/edward/bin/python
"""
This script makes fig 5
"""
from shared_functions import *

# ##### Figure 5 #####
def term_2(_df, lai, vpd):
  """calculates term 2"""
  atmos = {'gamma' : _df.gamma.mean(), 'c_a' : _df.c_a.mean(),\
           'vpd' : vpd}
  if _df.uwue.std() > 1.e-8:
    print('error, uWUE is variable: %f!!!!' % _df.uwue.std())
  elif _df.g1.std() > 1.e-8:
    print('error, g1 is variabile: %f!!!!!' % _df.g1.std())
  canopy = {'uwue' : _df.uwue.mean(), 'g1' : _df.g1.mean()}
  return pm.CP/_df.r_moist.mean() +calc.leaf_vpd(atmos, canopy, lai)


def term_2_approx(_df, lai, vpd, order=4):
  """calculates term 2"""
  atmos = {'gamma' : _df.gamma.mean(), 'c_a' : _df.c_a.mean(),\
           'vpd' : vpd}
  if _df.uwue.std() > 1.e-8:
    print('error, uWUE is variable: %f!!!!' % _df.uwue.std())
  elif _df.g1.std() > 1.e-8:
    print('error, g1 is variabile: %f!!!!!' % _df.g1.std())
  canopy = {'uwue' : _df.uwue.mean(), 'g1' : _df.g1.mean()}
  if order == 4:
    return pm.CP/_df.r_moist.mean() \
      -atmos['gamma']*atmos['c_a']*pm.LV/\
      (lai*1.6*pm.R_STAR*canopy['uwue'])\
      *(1./canopy['g1'] - 3.*np.sqrt(atmos['vpd'])/(2.*canopy['g1']**2)
        + 2.*atmos['vpd']/canopy['g1']**3\
        - 5.*np.sqrt(atmos['vpd'])**3/(2.*canopy['g1']**4))
  elif order == 2:
    return pm.CP/_df.r_moist.mean() \
      -atmos['gamma']*atmos['c_a']*pm.LV/\
      (lai*1.6*pm.R_STAR*canopy['uwue'])\
      *(1./canopy['g1'] - 3.*np.sqrt(atmos['vpd'])/(2.*canopy['g1']**2))
  else:
    print('error uncoded order number %d' % order)
    return


grouped = df.groupby('pft')
print('mean lai: %5.2f' % grouped.lai.mean().mean())
print('mean vpd: %5.2f' % grouped.vpd.mean().mean())

for key in ['lai', 'vpd', 'g1', 'uwue_norm']:
  print('cv %s: %5.2f' % (key,\
                          grouped[key].mean().std()/grouped[key].mean().mean()))

def get_pft(_df):
  return _df['pft'].iloc[0]

names = grouped.apply(get_pft)
plt.figure()
plt.plot(grouped.g1.mean(), grouped.uwue.mean(), 'k*')
for name, x, y in zip(names, grouped.g1.mean(), grouped.uwue.mean()):
  plt.annotate(name, xy=(x, y))
plt.savefig('%s/temp/garb.png' % os.environ['PLOTS'])

#now for figure 6 split product into two

### FIGURE 5 ###
# look at what controls variability between pfts
def first_half(_df, lai):
  """calcs first half of term 3"""
  return -_df.gamma.mean()*df.c_a.mean()/\
(lai*1.6*pm.R_STAR*_df.uwue_norm.mean())

def second_half(_df, vpd):
  """calcs second half of term 3"""
  _g1 = _df.g1.iloc[0]
  return -(2.*_g1 + np.sqrt(vpd))/\
    (2.*(_g1 + np.sqrt(vpd))**2)

def et_min_vpd(_df, lai):
  """calculates theoretical max vpd as functoin of -df and lai"""
  c3 = pm.CP/_df.r_moist
  c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
  c2 = _df.g1
  sqrt_vpd = (c1 + np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3)
  try:
    sqrt_vpd[sqrt_vpd < 0.] = np.nan
  except TypeError:
    if sqrt_vpd < 0.:
      sqrt_vpd = np.nan
  return sqrt_vpd**2

# def et_min_vpd1(_df, lai):
#   """calculates theoretical max vpd as functoin of -df and lai"""
#   """note below is only valid for negative x, which we don't have"""
#   c3 = pm.CP/_df.r_moist
#   c1 = _df.gamma*_df.c_a/(lai*pm.R_STAR*1.6*_df.uwue_norm)
#   c2 = _df.g1
#   return ((c1 - np.sqrt(c1 + 8.*c2*c3)*np.sqrt(c1)-4.*c2*c3)/(4.*c3))**2


I = 5
def pft_leaf(_df, axs):
  """takes df and plots both halves of product in term 2"""
  global I
  vpd = np.linspace(_df.vpd.quantile(q=0.05), _df.vpd.quantile(q=0.95))
  lai = _df.lai.mean()
  p = axs[0].plot(vpd, term_2(_df, lai, vpd),\
                  label="%s: $\overline{LAI}$=%4.2f, \n uWUE=%4.2f, g1=%4.1f"\
                  % (_df.pft.iloc[0],\
                     lai, _df.uwue_norm.iloc[0],  _df.g1.iloc[0]))
  axs[0].plot(vpd, term_2_approx(_df, lai, vpd, order=4), linestyle='dashed',\
              color=p[0].get_color())
  #print('axlim: ',axs[0].get_ylim())
  # axs[3].plot(vpd, second_half(_df, vpd),\
  #             label='%s, g1 = %4.1f' % (_df.pft.iloc[0], _df.g1.iloc[0]))
  ptiles = np.array([_df.vpd.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[1].plot(vpd, np.ones(vpd.shape)*I)
  axs[1].plot(ptiles, np.ones(ptiles.shape)*I, 'k*')
  axs[-1].plot(np.ones(vpd.shape)*I, vpd)
  axs[-1].plot(np.ones(ptiles.shape)*I, ptiles, 'k*')

  # # axs[2].plot(ptiles, term_2(_df, lai, ptiles), 'k*')
  # axs[3].plot(ptiles, second_half(_df, ptiles), 'k*')

  lai = np.linspace(_df.lai.quantile(q=0.05), _df.lai.quantile(q=0.95))
  vpd = _df.vpd.mean()
  _mean_df = _df.mean()
  et_min = et_min_vpd(_mean_df, lai)
  # et_min[et_min > _df.vpd.quantile(0.99)] = np.nan
  # et_min[et_min > 4001.] = np.nan
  axs[2].plot(lai, et_min,\
              label=r"PFT = %s, uWUE=%4.2f, g1=%4.1f"\
              % (_df.pft.iloc[0],\
                 _df.uwue_norm.iloc[0], _df.g1.iloc[0]))
  ptiles = np.array([_df.lai.quantile(q=_p/100.)\
                     for _p in [25., 50., 75.]])
  axs[3].plot(lai, np.ones(lai.shape)*I)
  axs[3].plot(ptiles, np.ones(ptiles.shape)*I, 'k*')
  # axs[1].plot(lai, first_half(_df, lai),\
  #             label='%s, uWUE = %4.2f'\
  #             % (_df.pft.iloc[0], _df.uwue_norm.iloc[0]))
  # # axs[0].plot(ptiles, term_2(_df, ptiles, vpd), 'k*')
  # axs[1].plot(ptiles, first_half(_df, ptiles), 'k*')
  # now second half
  I -= 1
  return

fig = plt.figure()
fig.set_figheight(fig.get_figheight()*2.5)
#fig.set_figwidth(fig.get_figwidth()*2)
#axs = [fig.add_subplot(2, 2, i+1) for i in range(4)]
_axs = [fig.add_subplot(2, 1, i+1) for i in range(2)]
axs = []
for _ax in _axs:
  divider = make_axes_locatable(_ax)
  axs.append(_ax)
  axs.append(divider.append_axes("bottom", size="20%", pad=0.0, sharex=_ax))
axs.append(divider.append_axes("left", size="20%", pad=0.0, sharey=axs[2]))
I = 5

for pft in ['CRO', 'DBF', 'GRA', 'ENF', 'CSH']:
  _df = df.loc[df.pft == pft, :]
  pft_leaf(_df, axs)

# df.groupby('pft').apply(pft_leaf, axs)
axs[1].set_xlabel('VPD (Pa)')
axs[3].set_xlabel('LAI')
axs[0].set_ylabel(paren_string)
plt.setp(axs[2].get_yticklabels(), visible=False)
axs[-1].set_ylabel(r'VPD')#$_{ETmin}$')
axs[2].set_title('VPD where ET = Min(ET) '\
                 r'($\frac{\partial \; ET}{\partial \; D} = 0$)')
axs[0].plot(axs[2].get_xlim(), [0., 0.], 'k--', linewidth=0.2)
# axs[2].set_ylim([0., np.around(df.vpd.quantile(q=0.95), decimals=-2)])
axs[2].set_ylim([0., 4000.])
axs[1].set_ylim([0.5,5.5])
axs[3].set_ylim([0.5,5.5])
axs[1].get_yaxis().set_visible(False)
axs[3].get_yaxis().set_visible(False)
axs[-1].set_xlim([0.5,5.5])
axs[-1].get_xaxis().set_visible(False)

axs[2].text(0.35, 250., r'$\frac{\partial \; ET}{\partial \; D} < 0$',\
            fontsize=15)
axs[2].text(1.4, 3500., r'$\frac{\partial \; ET}{\partial \; D} > 0$',\
            fontsize=15)
# axs[3].set_xlim(axs[2].get_xlim())
# axs[1].set_xlim(axs[0].get_xlim())

# axs[1].set_ylabel(r'-$\frac{\gamma c_s }{LAI \; 1.6 \; R\;  uWUE }$')
# axs[2].set_xlabel('VPD (Pa)')
# axs[3].set_xlabel('VPD (Pa)')
# axs[2].set_ylabel(paren_string)
# axs[3].set_ylabel(r'-$\frac{2 g_1 + \sqrt{D}}{2 (g_1 + \sqrt{D})^2}$')
axs[0].set_ylim((-1.8058955891452384, 0.87408219563044132))
for ax in axs[:1]:
  h, l = ax.get_legend_handles_labels()
  ax.legend(h, l, loc='best', fontsize=9)
# plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../doc/paper/fig05.pdf')
