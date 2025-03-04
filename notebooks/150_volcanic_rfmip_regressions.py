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

# %%
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
import sys
import pandas as pd
import os
#import urllib
from scipy.stats import linregress
import requests

# %%
pl.rcParams['figure.figsize'] = (16/2.54, 16/2.54)
pl.rcParams['font.size'] = 11
pl.rcParams['font.family'] = 'Arial'
pl.rcParams['xtick.direction'] = 'out'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['figure.dpi'] = 96


# %%
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
check_and_download(
    '../data_input/input4mips/solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc',
    'http://aims3.llnl.gov/thredds/dodsC/user_pub_work/input4MIPs/CMIP6/CMIP/SOLARIS-HEPPA/SOLARIS-HEPPA-3-2/atmos/mon/multiple/gm/v20170103/solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc'
)

# %%
nc_future = Dataset('../data_input/input4mips/solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc')
tsi_future = nc_future.variables['tsi'][:(251*12)]
nc_future.close()

# %%
years = np.arange(1850, 2101, dtype=int)
isleap = np.zeros(251)
isleap[np.logical_and(years%4==0, np.logical_or(years%100!=0, years%400==0))] = 1

# %%
tsi = np.zeros(251)
for i, year in enumerate(years):
    weights = [31,28+isleap[i],31,30,31,30,31,31,30,31,30,31]
    tsi[i] = np.average(tsi_future[(i)*12:(1+i)*12], weights=weights)

# %%
tsi_baseline = np.mean(tsi[:24])
tsi_baseline

# %%
pl.plot(np.arange(1850, 2101), tsi)
pl.axhline(tsi_baseline, ls=':', color='k')

# %%
solar_erf = 0.25 * (tsi-tsi_baseline) * 0.71 * 0.72  # include RA
pl.plot(np.arange(1850, 2101), solar_erf)
pl.axhline(0, ls=':', color='k')

# %%
# pre-processed RFMIP-histnat runs. Method pretty similar to notebook that does the ERF for aerosols
rfmip_tier2 = pd.read_csv('../data_input/RFMIP-ERF-tier2.csv', index_col=0)
natural = rfmip_tier2[['CanESM5 NAT', 'GFDL-CM4 NAT', 'GISS-E2-1-G NAT', 'HadGEM3-GC31-LL NAT', 'IPSL-CM6A-LR NAT', 'MIROC6 NAT', 'NorESM2-LM NAT']]
natural

# %%
natural.loc[1850:2015].mean().mean()

# %%
np.mean(solar_erf[:165])

# %%
pl.plot(natural)
pl.axhline(0, ls=':', color='k')

# %%
pl.plot(natural.mean(axis=1) - solar_erf)
pl.axhline(0, ls=':', color='k')

# %%
volc = natural - solar_erf[:,None]
volc

# %%
pl.plot(volc)

# %%
volc.loc[1992,:]

# %%
# This dataset is processed from the Toohey and Sigl / CMIP6 / GloSSAC data - see supplement for info
# full crunching notebook will be released with AR6 :)
saod = pd.read_csv("../data_input/volcanic_sAOD_monthly_-50001-201912.csv", index_col=0)

# %%
saod_annual = np.zeros(165)
for year in range(1850,2015):
    saod_annual[year-1850] = saod.loc[year:year+1,:].mean()
saod_annual = pd.DataFrame(
    {
        'year' : np.arange(1850,2015,dtype=int),
        'sAOD' : saod_annual
    }
)
saod_annual.set_index('year',inplace=True)
saod_annual

# %%
saod_annual.mean()

# %%
volc_saod = volc.join(saod_annual-saod_annual.mean())

# %%
volc_saod

# %%
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'CanESM5 NAT'].values)

# %%
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'CanESM5 NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'CanESM5 NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'GFDL-CM4 NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'GFDL-CM4 NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'GISS-E2-1-G NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'GISS-E2-1-G NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'HadGEM3-GC31-LL NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'HadGEM3-GC31-LL NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'HadGEM3-GC31-LL NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'HadGEM3-GC31-LL NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'IPSL-CM6A-LR NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'IPSL-CM6A-LR NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'MIROC6 NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'MIROC6 NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
pl.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'NorESM2-LM NAT'].values)
sl, ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'NorESM2-LM NAT'].values)
pl.plot(np.linspace(0,0.12,100), sl*np.linspace(0,0.12,100)+ic)
sl

# %%
#pl.rcParams['font.size']=16
#fig, ax=pl.subplots(figsize=(9, 9))
fig, ax=pl.subplots(figsize=(16/2.54, 16/2.54))
slope = np.zeros(7)
colors = {
    'CanESM5'        : 'red',#'#1e4c24',
    'E3SM'           : 'darkorange',
    'GFDL-ESM4'      : 'yellowgreen', 
    'GFDL-CM4'       : 'yellow',#'green',
    'GISS-E2-1-G'    : 'green',#'#771d7b',
    'HadGEM3-GC31-LL': 'turquoise',
    'IPSL-CM6A-LR'   : 'teal',
    'MIROC6'         : 'blue',#b85fb7',
    'MRI-ESM2-0'     : 'blueviolet',
    'NorESM2-LM'     : 'purple',#'red',
    'UKESM1-0-LL'    : 'crimson',
}
slope = {}
for model in ['CanESM5', 'GFDL-CM4', 'GISS-E2-1-G', 'HadGEM3-GC31-LL', 'IPSL-CM6A-LR', 'MIROC6', 'NorESM2-LM']:
    ax.scatter(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'%s NAT' % model].values, color=colors[model], alpha=0.3)
    slope[model], ic, _,_,_ = linregress(volc_saod.loc[:2014,'sAOD'].values, volc_saod.loc[:2014,'%s NAT' % model].values)
    ax.plot(np.linspace(0,0.12,100), slope[model]*np.linspace(0,0.12,100)+ic, color=colors[model], label=r'%s $%4.1f \tau$' % (model, slope[model]))
ax.legend()
ax.set_ylabel('Volcanic effective radiative forcing, W m$^{-2}$')
ax.set_xlabel('Stratospheric aerosol optical depth anomaly w.r.t. 1850-2014 mean')
ax.set_title('Volcanic forcing in RFMIP piClim-histnat')
ax.axhline(0, color='k', ls=':')
ax.axvline(0, color='k', ls=':')
pl.tight_layout()
pl.savefig('../figures/figureS8.png', dpi=300)
pl.savefig('../figures/figureS8.pdf')

# %%
slope

# %%
dummy = np.zeros(7)
for i, value in enumerate(slope.values()):
    dummy[i]=value
    
dummy.mean()

# %%
