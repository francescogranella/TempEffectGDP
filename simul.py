import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D

import matplotlib

import context
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import signal
from tqdm import tqdm

context.pdsettings()

os.chdir(r'C:\Users\Granella\Dropbox (CMCC)\PhD\Research\TempEffectGDP')

# Global seed
np.random.seed(9874)

def filter(v, period):
    b, a = signal.butter(4, 1/period, output='ba', btype='lowpass', fs=1)
    return signal.filtfilt(b, a, v)

# Set up filter
periods = [0, 3,5,10,15]


# Read data
df = pd.read_csv('Data/temp_series.csv')
df.columns = [f's{x}' for x in np.arange(5000) + 1]

def gen_g(df):
    # Growth rate
    base_g = 0.01
    g_sd = 0.005
    beta = gamma = -0.05
    e = np.random.normal(0, g_sd, len(df))

    df['g_gamma'] = base_g + gamma * df.ts + e
    df['g_beta'] = base_g + beta * (df.ts - df.ts.shift(+1)) + e
    df.dropna(inplace=True)  # Drop first period
    return df

def gen_g_2(df):
    # Growth rate
    base_g = 0.01
    g_sd = 0.005
    beta = gamma = -0.05
    e = np.random.normal(0, g_sd, len(df))

    df['g_gamma'] = base_g + gamma * df.ts + 0.5*gamma*df.ts.shift(+1) #+ 0.5*gamma*df.ts.shift(+2) + e
    df['g_beta'] = base_g + beta * (df.ts - df.ts.shift(+1)) + e
    df.dropna(inplace=True)
    return df

def gen_g_3(df):
    # Growth rate
    base_g = 0.01
    g_sd = 0.005
    beta = gamma = -0.05
    e = np.random.normal(0, g_sd, len(df))

    df['g_gamma'] = base_g + gamma * df.ts - 0.5*gamma*df.ts.shift(+1) #- 0.5*gamma*df.ts.shift(+2) + e
    df['g_beta'] = base_g + beta * (df.ts - df.ts.shift(+1)) + e

    df.dropna(inplace=True)
    return df


def gen_g_e(df, f):
    # Growth rate
    base_g = 0.01
    g_sd = 0.005
    beta = gamma = -0.05
    e = np.random.normal(0, g_sd, len(df))
    e_ts = np.random.normal(0, np.abs(f*gamma), len(df))
    _ts = df.ts + e_ts
    df['g_gamma'] = base_g + gamma * _ts
    df['g_beta'] = base_g + beta * (_ts - _ts.shift(+1)) + e

    df.dropna(inplace=True)
    return df


def detrend(varname, _df):
    X = pd.DataFrame(data=sm.add_constant(np.c_[_df.index,_df.index**2]), index=_df.index)
    return sm.OLS(_df[varname], X, missing='drop').fit().resid


# For every period: pass filter, compute the ratio, estimate the OLS coefficients.
# Store ratio and coefficient.
def filter_reg(df):
    data = []
    for p in periods:
        if p==0:
            ts_filtered = df.ts
            ratio = 1
            X = sm.add_constant(pd.DataFrame(pd.concat([ts_filtered, ts_filtered.shift(+1)], axis=1)))
        else:
            ts_filtered = filter(df.ts, p)
            ratio = np.median(df.ts / ts_filtered)
            X = sm.add_constant(ts_filtered)
        gamma_hat = sm.OLS(df.g_gamma, X, missing='drop').fit().params[1]
        beta_hat = sm.OLS(df.g_beta, X, missing='drop').fit().params[1]
        data.append([beta_hat, gamma_hat, ratio])
    return data



data = []
data2 = []
data3 = []
for col in tqdm(df.columns[:100]):
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    data.append(filter_reg(_df))
    # g 2
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g_2(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    data2.append(filter_reg(_df))
    # g 3
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g_3(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    data3.append(filter_reg(_df))


def densityplot(df):
    plt.hist(df.gamma_hat, bins=50, histtype='step', label='gamma_hat', density=True)
    plt.hist(df.gamma_hat_adj, bins=50, histtype='step', label='gamma_hat_adj', density=True)
    plt.axvline(-1)
    plt.legend()
    plt.show()


fig, ax = plt.subplots()
ax.grid(zorder=0)
ax.axhline(0, color='black', zorder=0)
ax.axhline(-0.05, color='black', zorder=0)
for i,p in enumerate(periods):
    r = pd.DataFrame(np.array(data)[:,i,:], columns=['beta_hat', 'gamma_hat', 'ratio'])
    r['gamma_hat_adj'] = r.gamma_hat / r.ratio
    r['beta_hat_adj'] = r.beta_hat / r.ratio
    _m = r.mean()
    _sd = r.std()
    ax.errorbar(i-0.1, _m['gamma_hat'], yerr=1.96*_sd['gamma_hat'], fmt='.k', color=cm.tab10(0))
    ax.errorbar(i-0.1, _m['beta_hat'], yerr=1.96*_sd['beta_hat'], fmt='.k', color=cm.tab10(1))
    ax.errorbar(i, _m['gamma_hat_adj'], yerr=1.96*_sd['gamma_hat_adj'], fmt='.k', color=cm.tab10(2))
    ax.errorbar(i, _m['beta_hat_adj'], yerr=1.96*_sd['beta_hat_adj'], fmt='.k', color=cm.tab10(3))
custom_lines = [Line2D([0], [0], color=cm.tab10(0), lw=2),
                Line2D([0], [0], color=cm.tab10(1), lw=2),
                Line2D([0], [0], color=cm.tab10(2), lw=2),
                Line2D([0], [0], color=cm.tab10(3), lw=2)]
ax.legend(custom_lines, [r'$\hat{\gamma}$', r'$\hat{\beta}$', r'$\hat{\gamma}$ adj', r'$\hat{\beta}$ adj'],
          loc='upper center', bbox_to_anchor=(0.5,-0.1),
          ncol=2)
ax.set_xticklabels([0] + periods)
ax.xaxis.set_major_locator(plt.MaxNLocator(len(periods)))
plt.title(r'Adjustment of $\beta$ and $\gamma$ effects')
plt.tight_layout()
plt.savefig('img/Adjustment.png')
plt.show()

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
for j in range(3):
    _data = [data, data2, data3][j]
    ax = axs[j]
    ax.grid(zorder=0)
    ax.axhline(0, color='black', zorder=0)
    ax.axhline(-0.05, color='black', zorder=0)
    for i,p in enumerate(periods):
        r = pd.DataFrame(np.array(_data)[:,i,:], columns=['beta_hat', 'gamma_hat', 'ratio'])
        r['gamma_hat_adj'] = r.gamma_hat / r.ratio
        r['beta_hat_adj'] = r.beta_hat / r.ratio
        _m = r.mean()
        _sd = r.std()
        ax.errorbar(i-0.1, _m['gamma_hat'], yerr=1.96*_sd['gamma_hat'], fmt='.k', color=cm.tab10(0))
        ax.errorbar(i-0.1, _m['beta_hat'], yerr=1.96*_sd['beta_hat'], fmt='.k', color=cm.tab10(1))
        ax.errorbar(i, _m['gamma_hat_adj'], yerr=1.96*_sd['gamma_hat_adj'], fmt='.k', color=cm.tab10(2))
        ax.errorbar(i, _m['beta_hat_adj'], yerr=1.96*_sd['beta_hat_adj'], fmt='.k', color=cm.tab10(3))
    custom_lines = [Line2D([0], [0], color=cm.tab10(0), lw=2),
                    Line2D([0], [0], color=cm.tab10(1), lw=2),
                    Line2D([0], [0], color=cm.tab10(2), lw=2),
                    Line2D([0], [0], color=cm.tab10(3), lw=2)]
ax.legend(custom_lines, [r'$\hat{\gamma}$', r'$\hat{\beta}$', r'$\hat{\gamma}$ adj', r'$\hat{\beta}$ adj'])
ax.set_xticklabels([0] + periods)
ax.xaxis.set_major_locator(plt.MaxNLocator(len(periods)))
plt.savefig('img/Lags_simulation.png')
plt.show()


np.random.seed(1233)
data = []
data_25 = []
data_50 = []
data_100 = []
e_25 = []
e_50 = []
e_100 = []
for col in tqdm(df.columns[:100]):
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    data.append(filter_reg(_df))
    #
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    e_ts = np.random.normal(0, 0.25*_df.ts.std(), len(_df))
    _df.ts = _df.ts + e_ts
    data_25.append(filter_reg(_df))
    e_25.append(e_ts)
    #
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta', _df)
    e_ts = np.random.normal(0, 0.5*_df.ts.std(), len(_df))
    _df.ts = _df.ts + e_ts
    data_50.append(filter_reg(_df))
    e_50.append(e_ts)
    #
    _df = df[[col]]
    _df.columns = ['ts']
    _df = gen_g(_df)
    _df['ts'] = detrend('ts', _df)
    _df['g_gamma'] = detrend('g_gamma', _df)
    _df['g_beta'] = detrend('g_beta',+ _df)
    e_ts = np.random.normal(0, 1*_df.ts.std(), len(_df))
    _df.ts = _df.ts + e_ts
    data_100.append(filter_reg(_df))
    e_100.append(e_ts)

fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12))
for j in range(4):
    _data = [data, data_25, data_50, data_100][j]
    _e = [[0], e_25, e_50, e_100][j]
    ax = axs[j]
    ax.grid(zorder=0)
    ax.axhline(0, color='black', zorder=0)
    ax.axhline(-0.05, color='black', zorder=0)
    ax.set_title('SD of noise: ' + str(np.round(np.mean(np.abs(_e)), 3)) + '°, ' + ['0', '25%', '50%', '100%'][j] + ' of SD of detrended simulated temperature', fontsize=10)
    if j==3:
        ax.set_xticklabels([0] + periods)
    else:
        ax.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
    # ax.set_title('SD of noise: ' + ['0', '25%', '50%', '100%'][j] + ' of SD of detrended simulated temperature', fontsize=10)
    for i,p in enumerate(periods):
        r = pd.DataFrame(np.array(_data)[:,i,:], columns=['beta_hat', 'gamma_hat', 'ratio'])
        r['gamma_hat_adj'] = r.gamma_hat / r.ratio
        r['beta_hat_adj'] = r.beta_hat / r.ratio
        _m = r.mean()
        _sd = r.std()
        ax.errorbar(i-0.1, _m['gamma_hat'], yerr=1.96*_sd['gamma_hat'], fmt='.k', color=cm.tab10(0))
        ax.errorbar(i-0.1, _m['beta_hat'], yerr=1.96*_sd['beta_hat'], fmt='.k', color=cm.tab10(1))
        ax.errorbar(i, _m['gamma_hat_adj'], yerr=1.96*_sd['gamma_hat_adj'], fmt='.k', color=cm.tab10(2))
        ax.errorbar(i, _m['beta_hat_adj'], yerr=1.96*_sd['beta_hat_adj'], fmt='.k', color=cm.tab10(3))
        # ax.plot(i, _m['gamma_hat_adj'], color=cm.tab10(2))
    custom_lines = [Line2D([0], [0], color=cm.tab10(0), lw=2),
                    Line2D([0], [0], color=cm.tab10(1), lw=2),
                    Line2D([0], [0], color=cm.tab10(2), lw=2),
                    Line2D([0], [0], color=cm.tab10(3), lw=2)]
axs[3].legend(custom_lines, [r'$\hat{\gamma}$', r'$\hat{\beta}$', r'$\hat{\gamma}$ adj', r'$\hat{\beta}$ adj'],
              loc='upper center', bbox_to_anchor=(0.5,-0.15), ncol=2)
ax.xaxis.set_major_locator(plt.MaxNLocator(len(periods)))
# plt.title(r'$\beta$ and $\gamma$ effects: adjustment')
plt.tight_layout()
plt.savefig('img/Attenuation_bias.png')
plt.show()

#
# fig, ax = plt.subplots(nrows=4, figsize=(10,5), sharex=True, sharey=True)
# for i,p in enumerate(periods):
#     _df = pd.DataFrame(data=np.array(data)[:,i,:], columns=['beta_hat', 'gamma_hat', 'ratio', 'gamma_hat_adj'])
#     _df['beta_hat_adj'] = _df['beta_hat'] / _df['ratio']
#     ax[i].hist(_df.gamma_hat, bins=50, histtype='step', label='gamma_hat', density=True)
#     ax[i].hist(_df.gamma_hat_adj, bins=50, histtype='step', label='gamma_hat_adj', density=True)
#     ax[i].hist(_df.beta_hat, bins=50, histtype='step', label='beta_hat', density=True)
#     ax[i].hist(_df.beta_hat_adj, bins=50, histtype='step', label='beta_hat_adj', density=True)
#     ax[i].set_title(p)
#     ax[i].grid()
#     ax[i].axvline(-0.05, color='black')
# plt.subplots_adjust(wspace=0, hspace=0)
# handles, labels = ax[i].get_legend_handles_labels()
# fig.legend(handles, labels)
# plt.show()


def _gen_g(_df):
    beta = gamma = -1
    _df['g_gamma'] = gamma * _df.ts
    _df['g_gamma_2'] = gamma * _df.ts + 0.5*gamma*_df.ts.shift(+1) # - 0.5*gamma*_df.ts.shift(+2)
    _df['g_gamma_3'] = gamma * _df.ts - 0.5*gamma*_df.ts.shift(+1) # + 0.5*gamma*_df.ts.shift(+2) # + 0.5*gamma*_df.ts.shift(+3) + 0.5*gamma*_df.ts.shift(+4)
    _df['g_gamma_5'] = _df.g_gamma[1:].cumsum() / 4
    _df['g_beta'] = beta * (_df.ts - _df.ts.shift(+1))
    _df['g_beta_2'] = beta * (_df.ts - _df.ts.shift(+1)) + 0.5*beta*(_df.ts.shift(+1) - _df.ts.shift(+2)) # + 0.5*beta*(_df.ts.shift(+2) - _df.ts.shift(+3)) + 0.5*beta*(_df.ts.shift(+3) - _df.ts.shift(+4))
    _df['g_beta_3'] = beta * (_df.ts - _df.ts.shift(+1)) - 0.5*beta*abs(_df.ts.shift(+1)==1)
    _df['g_beta_4'] = beta * (_df.ts - _df.ts.shift(+1)) - 0.5*beta*(_df.ts.shift(+1) - _df.ts.shift(+2)) # + 0.5*beta*(_df.ts.shift(+2) - _df.ts.shift(+3))  + 0.5*beta*(_df.ts.shift(+3) - _df.ts.shift(+4))
    _df.dropna(inplace=True)
    return _df

ts = [0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0]
dummy_df = pd.DataFrame({'ts': ts})
dummy_df = _gen_g(dummy_df).reset_index(drop=True)
fig, axs = plt.subplots(nrows=3, figsize=(10,10))
axs[0].plot(dummy_df.index, dummy_df.ts, color='gray', label='Temperature $T$')
for i in [1, 2]:
    axs[i].plot(dummy_df.index, dummy_df.g_beta, color=cm.tab20c(0), label=r'Beta: $\qquad\qquad\qquad\qquad\qquad\qquad g_t = -\Delta T_t$')
    axs[i].plot(dummy_df.index, dummy_df.g_beta_2, color=cm.tab20c(1), label=r'Beta Lag (same sign): $\qquad\qquad g_t = -\Delta T_t -0.5\Delta T_{t-1}$', linestyle='dashed')
    axs[i].plot(dummy_df.index, dummy_df.g_beta_4, color=cm.tab20c(2), label=r'Beta Lag (opposite sign): $\qquad g_t = -\Delta T_t +0.5\Delta T_{t-1}$', linestyle='dotted')
    axs[i].plot(dummy_df.index, dummy_df.g_beta_3, color=cm.tab20c(13), label=r'Beta compound: $\qquad\qquad\qquad g_t = -\Delta T_t - 0.5I[T_{t-1}\neq0]$', linestyle='dashdot')
    axs[i].plot(dummy_df.index, dummy_df.g_gamma, color=cm.tab20c(4), label=r'Gamma: $\qquad\qquad\qquad\qquad\qquad g_t = -T_t$')
    axs[i].plot(dummy_df.index, dummy_df.g_gamma_2, color=cm.tab20c(5), label=r'Gamma Lag (same sign): $\qquad g_t = -T_t - 0.5T_{t-1}$', linestyle='dashed')
    axs[i].plot(dummy_df.index, dummy_df.g_gamma_3, color=cm.tab20c(6), label=r'Gamma Lag (opposite sign): $ g_t = -T_t + 0.5T_{t-1}$', linestyle='dotted')
    box = axs[0].get_position()
    # axs[1].legend(loc='upper left', bbox_to_anchor=(0,-0.15), ncol=2)
    dummy_df = dummy_df.apply(lambda x: x.cumsum())
# axs[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
axs[0].set_ylabel('Temperature $T$')
axs[1].set_ylabel('Growth')
axs[2].set_ylabel('Log( GDP )')
axs[2].legend(loc='upper left', bbox_to_anchor=(0,-0.15), ncol=2)
plt.tight_layout()
plt.savefig('img/Fig1b_extended.png')
plt.show()


def _gen_g_2(_df):
    beta = gamma = -1
    _df['g_gamma'] = gamma * _df.ts
    _df['g_beta'] = beta * (_df.ts - _df.ts.shift(+1))
    _df.dropna(inplace=True)
    return _df


# Fig 1a
df = pd.read_csv(r"C:\Users\Granella\Box\Long Run GDP Growth\TempEffectGDP\Data\Barro_UDel.csv")
df = df.loc[df.countrycode=='USA', ['year', 'UDel_pop_temp']]
df = df.dropna(how='any')
df = df[df.year>=1960]
df.columns = ['year', 'temp']
df['temp'] = sm.OLS(df.temp, sm.add_constant(df.year)).fit().resid
df = df.set_index('year')

fig, ax = plt.subplots(figsize=(10,5))
ax.set(yticklabels=[])  # remove the tick labels
ax.tick_params(left=False)  # remove the ticks
ax.plot(df.temp, label='Unfiltered')
ax.axhline(df.temp.mean(), xmin=0.1, xmax=0.95, color=cm.tab10(0))
for i, p in enumerate([3,5,10,15]):
    df[f'f{p}'] = filter(df.temp,p) - 0.5*(i+1)
    ax.plot(df[f'f{p}'], label=f'2-{p} years')
    ax.axhline(df[f'f{p}'].mean(), xmin=0.1, xmax=0.95,  color=cm.tab10(i+1))
ax.errorbar(df.index.min()-4, df.temp.mean(), yerr=0.5, color='black')
ax.text(df.index.min()-3.5,  df.temp.mean(), "1C", fontsize=12)
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), title='Filtered oscillations (periods)')
plt.legend(loc='center left', bbox_to_anchor=(1,0.5), title='Filtered oscillations (periods)', frameon=False)
plt.title('Temperature fluctutations, United States')
plt.tight_layout()
plt.savefig('img/Fig1a.png')
plt.show()

# Lanczos filter
def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    Args:
    window: int
        The length of the filter window.
    cutoff: float
        The cutoff frequency in inverse time steps.
    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

fig, ax = plt.subplots(figsize=(10,5))
ax.set(yticklabels=[])  # remove the tick labels
ax.tick_params(left=False)  # remove the ticks
ax.plot(df.temp, label='Unfiltered')
ax.axhline(df.temp.mean(), xmin=0.1, xmax=0.95, color=cm.tab10(0))
window = len(df)
for i, p in enumerate([3,5,10,15]):
    fw = low_pass_weights(window, 1. / p)
    df[f'f{p}'] = np.convolve(fw, df.temp, mode='same') - 0.5*(i+1)
    ax.plot(df[f'f{p}'], label=f'2-{p} years')
    ax.axhline(df[f'f{p}'].mean(), xmin=0.1, xmax=0.95,  color=cm.tab10(i+1))
ax.errorbar(df.index.min()-4, df.temp.mean(), yerr=0.5, color='black')
ax.text(df.index.min()-3.5,  df.temp.mean(), "1C", fontsize=12)
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), title='Filtered oscillations (periods)')
plt.legend(loc='center left', bbox_to_anchor=(1,0.5), title='Filtered oscillations (periods)', frameon=False)
plt.title('Temperature fluctutations, United States. Lanczos filter')
plt.tight_layout()
plt.savefig('img/Fig1a-lanczos.png')
plt.show()

# Fig 1b
import matplotlib
tss = [ [0,0,1,0,0,0,0,0],
        [0,0,1,1,0,0,0,0],
        [0,0,1,1,1,0,0,0],
        ]
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(10,7), sharex=True)
for c in range(3):
    dummy_df = pd.DataFrame({'ts': tss[c]})
    dummy_df = _gen_g_2(dummy_df).reset_index(drop=True)
    axs[0,c].plot(dummy_df.index, dummy_df.ts, color='gray', label='Temperature $T$')
    axs[1,c].plot(dummy_df.index, dummy_df.g_beta, color=cm.tab10(0), linestyle='dashed', label=r'Level effect')
    axs[1,c].plot(dummy_df.index, dummy_df.g_gamma, color=cm.tab10(1), linestyle='dotted', label=r'Growth effect')
    dummy_df = dummy_df.apply(lambda x: x.cumsum())
    axs[2,c].plot(dummy_df.index, dummy_df.g_beta, color=cm.tab10(0), linestyle='dashed', label=r'Level effect')
    axs[2,c].plot(dummy_df.index, dummy_df.g_gamma, color=cm.tab10(1), linestyle='dotted', label=r'Growth effect')
    p = [(4,0.05), (4,-0.20), (4,-0.2)]
    p = [.1, .4,.85]
    for r in range(3):
        yl, yu = axs[r,c].get_ylim()
        axs[r,c].axhline(0, color='black', zorder=0)
        axs[r,c].text(4.5, yl + p[r]*(yu-yl), 'baseline')
        axs[r,c].yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
        axs[r,c].xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
        axs[r,c].xaxis.grid(alpha=0.5)
# axs[1,1].legend(loc='lower left', bbox_to_anchor=(0,-0.15), ncol=2)
axs[1,1].legend(loc='upper center', bbox_to_anchor=(0.5,1.2), ncol=2)
axs[0,0].set_ylabel('Temperature')
axs[1,0].set_ylabel('GDP Growth')
axs[2,0].set_ylabel('Log( GDP )')
fig.supxlabel('Time')
plt.tight_layout()
plt.savefig('img/Fig1b.png')
plt.show()