# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Historical projections with the new aerosol coefficients and emissions and constraints updated to 2022 and the dust forcing

# %%
# todo: most of these imports are redundant!

import numpy as np
import scipy.stats as st
import pandas as pd
import matplotlib.pyplot as pl
import os
import urllib
import json
import wquantiles
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
from scipy.optimize import root
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from tqdm.auto import tqdm
from scipy.signal import savgol_filter
from netCDF4 import Dataset
import matplotlib.gridspec as gridspec
import random
import h5py
from zipfile import ZipFile
from climateforcing.twolayermodel import TwoLayerModel
import scipy


# %%
# hdf5 utilities
def save_dict_to_hdf5(dic, filename):
    """
    ....
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans


# get my data
def check_and_download(filepath, url):
    """Checks prescence of a file and downloads if not present.

    Inputs
    ------
        filepath : str
            filename to download to
        url :
            url to download from
    """
    if not os.path.isfile(filepath):
        urllib.request.urlretrieve(url, filepath)
    return


# %%
pl.rcParams['figure.figsize'] = (12/2.54, 12/2.54)
pl.rcParams['font.size'] = 8
pl.rcParams['font.family'] = 'Arial'
pl.rcParams['xtick.direction'] = 'out'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['figure.dpi'] = 96

# %% [markdown]
# ## IGCC ocean heat uptake
#
# the citation is https://www.earth-syst-sci-data-discuss.net/essd-2019-255/

# %%
df_eei = pd.read_csv('../data_input/IGCC2022_earth_energy_imbalance.csv', index_col=0)
OHCobs = (df_eei.loc[:2020.5, 'Total'] - df_eei.loc[1971.5, 'Total']).values
print(OHCobs[-1])
#OHCobs_u = np.sqrt(ohctopu**2 + ohcbotu**2 + atmoshu**2 + cryoshu**2 + landhcu**2)
#pl.fill_between(np.arange(1960.5, 2019), OHCobs - 2 * OHCobs_u, OHCobs + 2 * OHCobs_u, alpha=0.3)
#pl.plot(np.arange(1960.5,2019), OHCobs)

# %% [markdown]
# ## Load in the dust forcing

# %%
dust_nc = Dataset('../data_input/Dust_radiative_forcing_timeseries_N1e5.nc')
dust_forcing = np.zeros((273, 100000))
#dust_nc.variables['year'][:]
#dust_nc.variables['n'][:]
#dust_nc.variables['ensemble members'][:]
dust_forcing[100:, :] = dust_nc.variables['ensemble members'][:173, :]
dust_forcing[:101, :] = np.linspace(0, dust_forcing[100, :], 101, axis=0)
pl.plot(dust_forcing);
dust_nc.close()

# %% [markdown]
# ## Non-aerosol forcing is based on IGCC 2022

# %%
ssp245_allforcing = pd.read_csv('../data_input/IGCC2022_ERF_best_aggregates_1750-2022.csv')
baseline_forcing = ssp245_allforcing[:].copy()

baseline_forcing.drop(
    labels=['timebound_lower','timebound_upper'],
    axis='columns',
    inplace=True
)
baseline_forcing.set_index('time', inplace=True)
pd.set_option('display.max_rows', 999)
baseline_forcing

# %%
baseline_forcing.index

# %%
# Temperature (GMST) observations: IGCC 2022
temp = pd.read_csv('../data_input/IGCC2022_annual_averages.csv')['gmst'].values

# %%
# GSAT/GMST ratio. Following IPCC we assume this is 1.
blratio = 1
years  = np.arange(1850.5, 2023)
Tobs = blratio * temp
pl.plot(years, Tobs)
#pl.plot(np.arange(1750,1901), best_land)
print(np.mean(Tobs[:51]))
print(np.mean(Tobs[160:170]))
#print(blratio)
len(Tobs)


# %% [markdown]
# ## Simple experiment with ECS=3.7 and default Geoffroy params

# %%
def rmse(obs, mod):
    return np.sqrt(np.sum((obs-mod)**2)/len(obs))


# %%
# load in Geoffroy two layer model parameters (pre-calculated by Glen Harris)
params = pd.read_fwf('../data_input/scmpy2L_calib_n=44_eps=fit_v20200702.txt', sep=' ')
params.set_index('Model', inplace=True)
cmip6_models = list(params.index)
params.rename(columns={"F4x":'q4x', "Lambda":'lamg', "Cmix":'cmix', "Cdeep":'cdeep', "Gamma":'gamma_2l', "Epsilon":'eff'}, inplace=True)

# %%
#cmip6_models = list(params['gamma_2l']['model_data']['EBM-epsilon'].keys())
fig, ax = pl.subplots(2,3, figsize=(16/2.54, 12/2.54))

ax[0,0].hist(params['q4x'], bins=np.arange(5,11,0.5), density=True)
f_q4x = st.gaussian_kde(params['q4x'], bw_method='silverman')
ax[0,0].plot(np.linspace(5,11), f_q4x(np.linspace(5,11)), color='k', label='Kernel density')
ax[0,0].set_yticks([])
ax[0,0].set_xlabel('W m$^{-2}$')
ax[0,0].set_title(r'ERF from $F_{4}\times$CO$_2$ ($F_{4\times}$)', fontsize=8)
ax[0,0].set_xlim(5,11)

ax[0,1].hist(-params['lamg'], bins=np.arange(-2,0.2,0.2), density=True)
f_lamg = st.gaussian_kde(-params['lamg'], bw_method='silverman')
ax[0,1].plot(np.linspace(-2.4,0), f_lamg(np.linspace(-2.4,0)), color='k', label='Kernel density')
ax[0,1].set_yticks([])
ax[0,1].set_xlabel('W m$^{-2}$ K$^{-1}$')
ax[0,1].set_title(r'Climate feedback parameter ($\lambda$)', fontsize=8)
ax[0,1].set_xlim(-2.4,0)

ax[0,2].hist(params['eff'], bins=np.arange(0.4,2.2,0.2), density=True)
f_eff = st.gaussian_kde(params['eff'], bw_method='silverman')
ax[0,2].plot(np.linspace(0.4,2.2), f_eff(np.linspace(0.4,2.2)), color='k', label='Kernel density')
ax[0,2].set_yticks([])
ax[0,2].set_title(r'Efficacy of ocean heat uptake ($\epsilon$)', fontsize=8)
ax[0,2].set_xlim(0.4,2.2)

ax[1,0].hist(params['gamma_2l'], bins=np.arange(0.3,1.3,0.1), density=True)
f_gamma_2l = st.gaussian_kde(params['gamma_2l'], bw_method='silverman')
ax[1,0].plot(np.linspace(0.2,1.1), f_gamma_2l(np.linspace(0.2,1.1)), color='k', label='Kernel density')
ax[1,0].set_yticks([])
ax[1,0].set_xlabel('W m$^{-2}$ K$^{-1}$')
ax[1,0].set_title(r'Heat exchange ($\gamma$)', fontsize=8)
ax[1,0].set_xlim(0.2,1.1)

ax[1,1].hist(params['cmix'], bins=np.arange(5,11.5,0.5), density=True)
f_cmix = st.gaussian_kde(params['cmix'], bw_method='silverman')
ax[1,1].plot(np.linspace(5,11.5), f_cmix(np.linspace(5,11.5)), color='k', label='Kernel density')
ax[1,1].set_yticks([])
ax[1,1].set_xlabel('W yr m$^{-2}$ K$^{-1}$')
ax[1,1].set_title(r'Mixed-layer heat capacity (C)', fontsize=8)
ax[1,1].set_xlim(5,11.5)

ax[1,2].hist(params['cdeep'], bins=np.arange(0,420, 20), density=True, label='CMIP6 models')
f_cdeep = st.gaussian_kde(params['cdeep'], bw_method='silverman')
ax[1,2].plot(np.arange(0,400), f_cdeep(np.arange(0,400)), color='k', label='Kernel density')
ax[1,2].set_yticks([])
ax[1,2].set_xlabel('W yr m$^{-2}$ K$^{-1}$')
ax[1,2].set_title(r'Deep ocean heat capacity (C$_0$)', fontsize=8)
ax[1,2].set_xlim(0,400)
ax[1,2].legend()

pl.figtext(0.015,0.775,'Probability density', rotation=90, va='center', ha='center')
pl.figtext(0.015,0.275,'Probability density', rotation=90, va='center', ha='center')

fig.tight_layout(rect=[0.015,0,1,1])

# %%
# construct correlation matrix
pd.set_option('display.precision', 4)
params.corr()

# %%
pd.set_option('display.precision', 2)
geoff_ecs_data = np.zeros(((len(cmip6_models))))
geoff_ecs_data = params['q4x'].values/params['lamg'].values/2
geoff_df_display = params.copy()
geoff_df_display['ECS'] = geoff_ecs_data
geoff_df_display.sort_index()

# %%
samples = 100000

kde = st.gaussian_kde(params.T)
geoff_sample = kde.resample(size=int(samples*1.1), seed=3170812)
# remove unphysical combinations
geoff_sample[:,geoff_sample[0,:] <= 0] = np.nan
geoff_sample[:,geoff_sample[1,:] <= 0.2] = np.nan
geoff_sample[:,geoff_sample[2,:] <= 0] = np.nan
geoff_sample[:,geoff_sample[3,:] <= 0] = np.nan
geoff_sample[:,geoff_sample[4,:] <= 0] = np.nan
geoff_sample[:,geoff_sample[5,:] <= 0] = np.nan
#geoff_sample = geoff_sample[~np.isnan(geoff_sample)]
mask = np.all(np.isnan(geoff_sample), axis=0)
geoff_sample = geoff_sample[:,~mask]
geoff_sample_df=pd.DataFrame(
    data=geoff_sample[:,:samples].T, columns=['q4x','lamg','cmix','cdeep','gamma_2l','eff']
)
geoff_sample_df

# %%
# fractional uncertainties on ERF - based on the FAIR code
seed    = 36572 
zscore = st.norm.ppf(0.95)

# update these ranges for AR6
unc_ranges = np.array([
    0.12,      # CO2
    0.20,      # CH4: updated value from etminan 2016
    0.14,      # N2O
    0.19,      # other WMGHGS
    0.50,      # O3
    1.00,      # stratospheric WV from CH4
    0.00,      # contrails (non-symmetric)
    0.00,      # black carbon on snow (non-symmetric)
    0.50,      # land use change
    0.25,      # volcanic
    0.50,      # solar (amplitude)
])/(zscore)

def opt(x, q05_desired, q50_desired, q95_desired):
    "x is (a, loc, scale) in that order."
    q05, q50, q95 = scipy.stats.skewnorm.ppf(
        (0.05, 0.50, 0.95), x[0], loc=x[1], scale=x[2]
    )
    return (q05 - q05_desired, q50 - q50_desired, q95 - q95_desired)

scale = st.norm.rvs(size=(samples,11), loc=np.ones((samples,11)), scale=np.ones((samples, 11)) * unc_ranges[None,:], random_state=seed)

lapsi_params = scipy.optimize.root(opt, [1, 1, 1], args=(0, 1, 2.25)).x
contrails_params = scipy.optimize.root(opt, [1, 1, 1], args=(19 / 57, 1, 98 / 57)).x

# contrails 
scale[:, 6] = scipy.stats.skewnorm.rvs(
    contrails_params[0],
    loc=contrails_params[1],
    scale=contrails_params[2],
    size=samples,
    random_state=3701585,
)

# lapsi
scale[:, 7] = scipy.stats.skewnorm.rvs(
    lapsi_params[0],
    loc=lapsi_params[1],
    scale=lapsi_params[2],
    size=samples,
    random_state=3701584,
)

scale_df = pd.DataFrame(
    data = scale,
    columns = ['CO2','CH4','N2O','halogen','O3','H2O_stratospheric','contrails','BC_on_snow','land_use','volcanic','solar']
)
scale_df

# %%
pl.hist(scale[:,7])

# %%
trend_solar = st.norm.rvs(size=samples, loc=0, scale=0.1/zscore, random_state=138294)

# %% [markdown]
# ## Get SLCFs from IGCC 2022

# %%
emissions_ceds_update = pd.read_csv('../data_input/IGCC2024_slcf_emissions_1750-2024.csv', index_col=0)
emissions_ceds_update

# %%
#emissions = pd.read_csv('../output_data/historical_slcf_emissions.csv', index_col='year')
emissions = emissions_ceds_update.drop(['CO','NMVOC','NOx','NH3'], axis=1)
emissions


# %%
def ari_linear_nobase(x, a0, a1, a2):
    ari = x[0] * a0 + x[1] * a1 + x[2] * a2
    return ari

def aci_log_nobase(x, beta, n0, n1, n2):
    aci = beta * np.log(1 + x[0] * n0 + x[1] * n1 + x[2] * n2)
    return aci

df = pd.read_csv('../data_output/corrected_ERFari_samples.csv')
ari_coeffs = df.values

df = pd.read_csv('../data_output/corrected_ERFaci_samples.csv')
aci_coeffs = np.exp(df.values)


# %%
# Use Ringberg aerosol priors, from script provided to me
def uniform1684(a,b,seed,samples=2000):
    interval = (b-a)+((b-a)/(84-16)*32)
    lower = a-((b-a)/(84-16)*16)
    return st.uniform.rvs(lower, interval, size=samples, random_state=seed)

dtau    = uniform1684(0.02,0.04,123,samples=samples)
tau     = uniform1684(0.13,0.17,124,samples=samples)
S_tau   = uniform1684(-27,-20,125,samples=samples)
RFari_cloudy = uniform1684(-0.1,0.1,126,samples=samples)
dR_dRatm = uniform1684(-0.3,-0.1,127,samples=samples)
dRatm_dtau = uniform1684(17,35,128,samples=samples)
c_tau    = uniform1684(0.59,0.71,129,samples=samples)
c_N      = uniform1684(0.19,0.29,130,samples=samples)
c_L      = uniform1684(0.21,0.29,131,samples=samples)
c_C      = uniform1684(0.59,1.07,132,samples=samples)
beta_N_tau = uniform1684(0.3,0.8,133,samples=samples)
beta_L_N   = uniform1684(-0.36,-0.011,134,samples=samples)
beta_C_N   = uniform1684(0,0.1,135,samples=samples)
S_N = uniform1684(-27,-26,136,samples=samples)
S_L = uniform1684(-56,-54,137,samples=samples)
S_C = uniform1684(-153,-91,138,samples=samples)

rfari = dtau*S_tau*(1-c_tau)+RFari_cloudy
rfari_adj = dtau*dR_dRatm*dRatm_dtau

dlntau = dtau/tau
deltan = dlntau * beta_N_tau

rfaci = dlntau*beta_N_tau*S_N*c_N
erfaci_L = dlntau*beta_N_tau*beta_L_N*S_L*c_L
erfaci_C = dlntau*beta_N_tau*beta_C_N*S_C*c_C

ERFari_scale = rfari + rfari_adj
ERFaci_scale = rfaci + erfaci_L + erfaci_C

# %%
fig,ax = pl.subplots(1,2, figsize=(18/2.54,12/2.54))
ax[0].hist(ERFari_scale, bins=np.arange(-1.2,0.05,0.05), density=True);
ax[0].set_title('ERFari prior');
ax[1].hist(ERFaci_scale, bins=np.arange(-4,0.5,0.2), density=True);
ax[1].set_title('ERFaci prior');

# %%
# Define our dicts
ERFari = {}
ERFaci = {}
temp = {}
ks = {}
ohc = {}
hflux = {}

ks['temp'] = {}
ks['ohc'] = {}
ks['multi'] = {}

# %%
# # load dicts
# ERFari = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/ERFari.h5')
# ERFaci = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/ERFaci.h5')
# temp   = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/temp.h5')
# ks     = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/knutti_score.h5')
# ohc    = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/ohc.h5')
# hflux  = load_dict_from_hdf5('/nfs/a65/pmcjs/AR6_tuning/aerosols/hflux.h5')

# %%
intvar = np.loadtxt('../data_output/piControl/internal_variability_piControl.txt')


# %%
def knutti_score(obs, mod, sigma_D=None):
    """
    obs: observations data: array of size (nyears,)
    mod: model data: array of size (nyears, nsamples)
    """
    samples = mod.shape[1]
    rm_d = np.ones(samples) * np.nan
    for i in range(samples):
        rm_d[i] = rmse(obs, mod[:, i])
    if sigma_D==None:
        sigma_D = np.nanmin(rm_d)
    veracity = np.exp(-rm_d**2/sigma_D**2) 
    ks_raw = veracity
    ks_raw[np.isnan(ks_raw)] = 0
    ks = ks_raw/np.sum(ks_raw)
    return ks


# %%
def weighted_percentile(a, w, q):
    if isinstance(q, (list, tuple, np.ndarray)):
        result = []
        for iq in q:
            result.append(wquantiles.quantile(a, w, iq))
    else:
        result = wquantiles.quantile(a, w, q)
    return result


# %%
def simple_weight(obs, mod, sigma_D):
    veracity = np.exp(-(mod-obs)**2/sigma_D**2)
    similarity = 1 # I see no good reason to change this #
    ks_raw = veracity/similarity
    ks_raw[np.isnan(ks_raw)] = 0
    ks = ks_raw/np.sum(ks_raw)
    return ks


# %% [markdown]
# ## CEDS emissions

# %%
emissions.loc[2005:2015, 'OC']

# %%
ERFari['CMIP6-constrained'] = np.zeros((273,samples))
for i in tqdm(range(samples)):
    ts2010 = np.mean(
        ari_linear_nobase(
            [
                emissions.loc[2005:2015, 'SO2'], 
                emissions.loc[2005:2015, 'BC'], 
                emissions.loc[2005:2015, 'OC']
            ],
            ari_coeffs[i,0], ari_coeffs[i,1], ari_coeffs[i,2]
        )
    )

    ts1850 = np.mean(
        ari_linear_nobase(
            [
                emissions.loc[1850, 'SO2'], 
                emissions.loc[1850, 'BC'], 
                emissions.loc[1850, 'OC']
            ],
            ari_coeffs[i,0], ari_coeffs[i,1], ari_coeffs[i,2]
        )
    )

    ts1750 = np.mean(
        ari_linear_nobase(
            [
                emissions.loc[1750, 'SO2'], 
                emissions.loc[1750, 'BC'], 
                emissions.loc[1750, 'OC']
            ],
            ari_coeffs[i,0], ari_coeffs[i,1], ari_coeffs[i,2]
        )
    )
    
    ERFari['CMIP6-constrained'][:,i] = (
        ari_linear_nobase(
            [
                emissions.loc[:2022, 'SO2'], 
                emissions.loc[:2022, 'BC'], 
                emissions.loc[:2022, 'OC']
            ], ari_coeffs[i,0], ari_coeffs[i,1], ari_coeffs[i,2]
        ) - ts1750
    ) / (ts2010 - ts1850) * (ERFari_scale[i]) + dust_forcing[:, i]

# %%
ERFaci['CMIP6-constrained'] = np.zeros((273,samples))
for i in tqdm(range(samples)):
    ts2010 = np.mean(
        aci_log_nobase(
            [
                emissions.loc[2005:2015,'SO2'], 
                 emissions.loc[2005:2015,'BC'],
                 emissions.loc[2005:2015,'OC']
            ], 1, aci_coeffs[i,0], aci_coeffs[i,1], aci_coeffs[i,2]
        )
    )
    
    ts1850 = aci_log_nobase(
        [
            emissions.loc[1850,'SO2'],
            emissions.loc[1850,'BC'],
            emissions.loc[1850,'OC']
        ], 1, aci_coeffs[i,0], aci_coeffs[i,1], aci_coeffs[i,2]
    )

    ts1750 = aci_log_nobase(
        [
            emissions.loc[1750,'SO2'],
            emissions.loc[1750,'BC'],
            emissions.loc[1750,'OC']
        ], 1, aci_coeffs[i,0], aci_coeffs[i,1], aci_coeffs[i,2]
    )
    
    ERFaci['CMIP6-constrained'][:,i] = (
        aci_log_nobase(
            [
                emissions.loc[:2022, 'SO2'],
                emissions.loc[:2022, 'BC'],
                emissions.loc[:2022, 'OC']
            ], 1, aci_coeffs[i,0], aci_coeffs[i,1], aci_coeffs[i,2]
        ) - ts1750
    ) / (ts2010-ts1850) * (ERFaci_scale[i])

# %%
fig, ax = pl.subplots(1,3,figsize=(19/2.54, 9.5/2.54))
ax[0].fill_between(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained'], 5, axis=1), np.percentile(ERFari['CMIP6-constrained'], 95, axis=1), color='0.75', lw=0);
ax[0].fill_between(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained'], 16, axis=1), np.percentile(ERFari['CMIP6-constrained'], 84, axis=1), color='0.5', lw=0);
ax[0].plot(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained'], 50, axis=1), color='k', zorder=10)
ax[0].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,754], color='cyan', label='Parameter set #754')
ax[0].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,1076], color='magenta', label='Parameter set #1076')
ax[0].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,18010], color='lime', label='Parameter set #18010')
ax[0].legend()
ax[0].set_xlim(1800,2022)
ax[0].set_title('ERFari')
ax[0].set_ylabel('W m$^{-2}$')

ax[1].fill_between(np.arange(1750,2023), np.percentile(ERFaci['CMIP6-constrained'], 5, axis=1), np.percentile(ERFaci['CMIP6-constrained'], 95, axis=1), color='0.75', lw=0, label='5-95% range');
ax[1].fill_between(np.arange(1750,2023), np.percentile(ERFaci['CMIP6-constrained'], 16, axis=1), np.percentile(ERFaci['CMIP6-constrained'], 84, axis=1), color='0.5', lw=0, label='16-84% range');
ax[1].plot(np.arange(1750,2023), np.percentile(ERFaci['CMIP6-constrained'], 50, axis=1), color='k', label='median', zorder=10)
ax[1].plot(np.arange(1750,2023), ERFaci['CMIP6-constrained'][:,754], color='cyan')
ax[1].plot(np.arange(1750,2023), ERFaci['CMIP6-constrained'][:,1076], color='magenta')
ax[1].plot(np.arange(1750,2023), ERFaci['CMIP6-constrained'][:,18010], color='lime')
ax[1].legend()
ax[1].set_xlim(1800,2022)
ax[1].set_title('ERFaci')

ax[2].fill_between(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained']+ERFaci['CMIP6-constrained'], 5, axis=1), np.percentile(ERFari['CMIP6-constrained']+ERFaci['CMIP6-constrained'], 95, axis=1), color='0.75', lw=0);
ax[2].fill_between(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained']+ERFaci['CMIP6-constrained'], 16, axis=1), np.percentile(ERFari['CMIP6-constrained']+ERFaci['CMIP6-constrained'], 84, axis=1), color='0.5', lw=0);
ax[2].plot(np.arange(1750,2023), np.percentile(ERFari['CMIP6-constrained']+ERFaci['CMIP6-constrained'], 50, axis=1), color='k', zorder=10)
ax[2].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,754]+ERFaci['CMIP6-constrained'][:,754], color='cyan', label='Ensemble 754')
ax[2].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,1076]+ERFaci['CMIP6-constrained'][:,1076], color='magenta', label='Ensemble 1076')
ax[2].plot(np.arange(1750,2023), ERFari['CMIP6-constrained'][:,18010]+ERFaci['CMIP6-constrained'][:,18010], color='lime', label='Ensemble 18010')
ax[2].set_xlim(1800,2022)
ax[2].set_title('Aerosol ERF')

ax[0].set_ylim(-3,0.2)
ax[1].set_ylim(-3,0.2)
ax[2].set_ylim(-3,0.2)
ax[0].axhline(0, ls=':', lw=0.5, color='k')
ax[1].axhline(0, ls=':', lw=0.5, color='k')
ax[2].axhline(0, ls=':', lw=0.5, color='k')

fig.tight_layout()
# pl.savefig('../figures/figure5.png', dpi=300)
# pl.savefig('../figures/figure5.pdf')

# %%
in_forcing = baseline_forcing.copy()
in_forcing.drop(
    [
        'aerosol-radiation_interactions', 
        'aerosol-cloud_interactions',
        'aerosol',
        'anthro',
        'nonco2wmghg',
        'minor',
        'total'
    ], axis=1, inplace=True
)
in_forcing
in_forcing = in_forcing * scale_df.loc[0,:]
in_forcing['solar'] = in_forcing['solar'] + np.linspace(0, trend_solar[0], 273)
in_forcing['aerosol-radiation_interactions'] = ERFari['CMIP6-constrained'][:273,0]
in_forcing['aerosol-cloud_interactions'] = ERFaci['CMIP6-constrained'][:273,0]
in_forcing['total'] = in_forcing.sum(axis=1)

# %%
temp['CMIP6-constrained'] = np.zeros((273, samples))
ohc['CMIP6-constrained'] = np.zeros((273, samples))
hflux['CMIP6-constrained'] = np.zeros((273, samples))
for i in tqdm(range(samples)):
    in_forcing = baseline_forcing.copy()
    in_forcing.drop(
        [
            'aerosol-radiation_interactions', 
            'aerosol-cloud_interactions',
            'aerosol',
            'anthro',
            'nonco2wmghg',
            'minor',
            'total'
        ], axis=1, inplace=True
    )
    in_forcing = in_forcing * scale_df.loc[i,:]
    in_forcing['solar'] = in_forcing['solar'] + np.linspace(0, trend_solar[i], 273)
    in_forcing['aerosol-radiation_interactions'] = ERFari['CMIP6-constrained'][:273,i]
    in_forcing['aerosol-cloud_interactions'] = ERFaci['CMIP6-constrained'][:273,i]
    in_forcing['total'] = in_forcing.sum(axis=1)
    scm = TwoLayerModel(
        extforce=in_forcing['total'],
        exttime=in_forcing.index,
        tbeg=1750,
        tend=2023,
        q2x=geoff_sample_df.loc[i,'q4x']/2,
        lamg=geoff_sample_df.loc[i,'lamg'],
        t2x=None,
        eff=geoff_sample_df.loc[i,'eff'],
        cmix=geoff_sample_df.loc[i,'cmix'],
        cdeep=geoff_sample_df.loc[i,'cdeep'],
        gamma_2l=geoff_sample_df.loc[i,'gamma_2l'],
        outtime=np.arange(1750.5,2023),
        dt=1
    )
    out = scm.run()
    temp['CMIP6-constrained'][:,i] = out.tg - np.mean(out.tg[100:150])
    ohc['CMIP6-constrained'][:,i] = out.ohc
    hflux['CMIP6-constrained'][:,i] = out.hflux

# %%
pl.fill_between(np.arange(1750,2023), np.percentile(temp['CMIP6-constrained'], 5, axis=1), np.percentile(temp['CMIP6-constrained'], 95, axis=1))
pl.plot(np.arange(1750,2023), np.median(temp['CMIP6-constrained'], axis=1), color='k')

# %%
ks['temp']['CMIP6-constrained'] = knutti_score(Tobs, temp['CMIP6-constrained'][100:273, :] + intvar[100:273,:samples], sigma_D=0.12)
    # unchanged sigma_D for temperature; slightly unsatisfactory since it's not symmetric
ks['ohc']['CMIP6-constrained'] = simple_weight(465.3, 10*(ohc['CMIP6-constrained'][270,:]-ohc['CMIP6-constrained'][221,:]), sigma_D=66)
    # ohc sigma_D range comes from fair-calibrate 1.4.1 which is based on IGCC 2022
ks['multi']['CMIP6-constrained'] = (ks['temp']['CMIP6-constrained']*ks['ohc']['CMIP6-constrained'])/(np.sum(ks['temp']['CMIP6-constrained']*ks['ohc']['CMIP6-constrained']))

# %%
print(weighted_percentile(ERFari['CMIP6-constrained'][272,:]+ERFaci['CMIP6-constrained'][272,:], ks['temp']['CMIP6-constrained'][:], [.05,.16,.5,.84,.95]))
print(weighted_percentile(ERFari['CMIP6-constrained'][272,:]+ERFaci['CMIP6-constrained'][272,:], ks['ohc']['CMIP6-constrained'][:], [.05,.16,.5,.84,.95]))
print(weighted_percentile(ERFari['CMIP6-constrained'][272,:]+ERFaci['CMIP6-constrained'][272,:], ks['multi']['CMIP6-constrained'][:], [.05,.16,.5,.84,.95]))

# %%
os.makedirs('../data_output/results/dust/', exist_ok=True)

# %%
save_dict_to_hdf5(ERFari, '../data_output/results/dust/ERFari.h5')
save_dict_to_hdf5(ERFaci, '../data_output/results/dust/ERFaci.h5')
save_dict_to_hdf5(temp, '../data_output/results/dust/temp.h5')
save_dict_to_hdf5(ks, '../data_output/results/dust/knutti_score.h5')
save_dict_to_hdf5(ohc, '../data_output/results/dust/ohc.h5')
save_dict_to_hdf5(hflux, '../data_output/results/dust/hflux.h5')

# %%
# Throw TCR into the mix
samples = 100000
tcr = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'] + geoff_sample_df[:samples]['eff']*geoff_sample_df[:samples]['gamma_2l'])
ecs = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'])

# %%
expts = ['CMIP6-constrained']

# %%
pc = {}
for expt in tqdm(expts):
    pc[expt] = {}
    for constraint in ['temp', 'ohc', 'multi']:
        pc[expt][constraint] = {}
        pc[expt][constraint]['ECS'] = {}
        pc[expt][constraint]['TCR'] = {}
        for metric in ['GSAT','OHC','ERFari','ERFaci','ERFaer']:
            pc[expt][constraint][metric] = {}
            for perc in ['5','16','50','84','95']:
                pc[expt][constraint][metric][perc] = np.zeros(273)
        (
            pc[expt][constraint]['ECS']['5'],
            pc[expt][constraint]['ECS']['16'],
            pc[expt][constraint]['ECS']['50'],
            pc[expt][constraint]['ECS']['84'],
            pc[expt][constraint]['ECS']['95']
        ) = (
            weighted_percentile(ecs, ks[constraint][expt], [.05,.16,.5,.84,.95])
        )
        (
            pc[expt][constraint]['TCR']['5'],
            pc[expt][constraint]['TCR']['16'],
            pc[expt][constraint]['TCR']['50'],
            pc[expt][constraint]['TCR']['84'],
            pc[expt][constraint]['TCR']['95']
        ) = weighted_percentile(tcr, ks[constraint][expt], [.05,.16,.5,.84,.95])
        for year in range(273):
            (
                pc[expt][constraint]['GSAT']['5'][year],
                pc[expt][constraint]['GSAT']['16'][year],
                pc[expt][constraint]['GSAT']['50'][year],
                pc[expt][constraint]['GSAT']['84'][year],
                pc[expt][constraint]['GSAT']['95'][year] 
            ) = weighted_percentile(temp[expt][year,:] + intvar[year,:samples], ks[constraint][expt], [.05,.16,.5,.84,.95])
            (
                pc[expt][constraint]['OHC']['5'][year],
                pc[expt][constraint]['OHC']['16'][year],
                pc[expt][constraint]['OHC']['50'][year],
                pc[expt][constraint]['OHC']['84'][year],
                pc[expt][constraint]['OHC']['95'][year] 
            ) = weighted_percentile(ohc[expt][year,:], ks[constraint][expt], [.05,.16,.5,.84,.95])
            (
                pc[expt][constraint]['ERFari']['5'][year],
                pc[expt][constraint]['ERFari']['16'][year],
                pc[expt][constraint]['ERFari']['50'][year],
                pc[expt][constraint]['ERFari']['84'][year],
                pc[expt][constraint]['ERFari']['95'][year] 
            ) = weighted_percentile(ERFari[expt][year,:], ks[constraint][expt], [.05,.16,.5,.84,.95])
            (
                pc[expt][constraint]['ERFaci']['5'][year],
                pc[expt][constraint]['ERFaci']['16'][year],
                pc[expt][constraint]['ERFaci']['50'][year],
                pc[expt][constraint]['ERFaci']['84'][year],
                pc[expt][constraint]['ERFaci']['95'][year] 
            ) = weighted_percentile(ERFaci[expt][year,:], ks[constraint][expt], [.05,.16,.5,.84,.95])
            (
                pc[expt][constraint]['ERFaer']['5'][year],
                pc[expt][constraint]['ERFaer']['16'][year],
                pc[expt][constraint]['ERFaer']['50'][year],
                pc[expt][constraint]['ERFaer']['84'][year],
                pc[expt][constraint]['ERFaer']['95'][year] 
            ) = weighted_percentile(ERFari[expt][year,:]+ERFaci[expt][year,:], ks[constraint][expt], [.05,.16,.5,.84,.95])

# %%
save_dict_to_hdf5(pc, '../data_output/results/dust/pc.h5')

# %%
