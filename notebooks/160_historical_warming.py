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
import os
import numpy as np
import scipy.stats as st
import glob
import matplotlib.pyplot as pl
import warnings
import h5py
from tqdm import tqdm_notebook
from scipy.interpolate import interp1d
warnings.simplefilter('ignore')

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
models = []

for path in glob.glob(os.path.join(*os.path.normpath('../data_input/cmip6/*').split(os.sep))):
    models.append(os.path.split(path)[-1])

models

# %%
# models = [
# 'ACCESS-CM2',     'CanESM5-CanOE',   'CNRM-CM6-1',    'EC-Earth3-LR',     'GISS-E2-1-G',     'INM-CM5-0',       'NESM3',
# 'ACCESS-ESM1-5',  'CAS-ESM2-0',      'CNRM-CM6-1-HR', 'EC-Earth3-Veg',    'GISS-E2-1-G-CC',  'IPSL-CM6A-LR',    'NorCPM1',
# 'AWI-CM-1-1-MR',  'CESM2',           'CNRM-ESM2-1',   'EC-Earth3-Veg-LR', 'GISS-E2-1-H',     'MIROC6',          'NorESM1-F',
# 'AWI-ESM-1-1-LR', 'CESM2-FV2',                        'FGOALS-f3-L',      'GISS-E2-2-G',     'MIROC-ES2L',      'NorESM2-LM',
# 'BCC-CSM2-MR',    'CESM2-WACCM',     'E3SM-1-0',      'FGOALS-g3',        'HadGEM3-GC31-LL', 'MPI-ESM-1-2-HAM', 'NorESM2-MM',
# 'BCC-ESM1',       'CESM2-WACCM-FV2', 'E3SM-1-1',      'FIO-ESM-2-0',      'HadGEM3-GC31-MM', 'MPI-ESM1-2-HR',   'SAM0-UNICON',
# 'CAMS-CSM1-0',    'CIESM',           'E3SM-1-1-ECA',  'GFDL-CM4',         'IITM-ESM',        'MPI-ESM1-2-LR',   'TaiESM1',
# 'CanESM5',        'CMCC-CM2-SR5',    'EC-Earth3',     'GFDL-ESM4',        'INM-CM4-8',       'MRI-ESM2-0',      'UKESM1-0-LL'
# ]

# %%
historical = {}
accepted_models = []
nyears = {}

for model in models:
    historical[model] = {}
    path_hist_tas  = glob.glob(
        os.path.join(
            '..', 'data_input', 'cmip6', model, 'historical', '*', 'tas.txt'
        )
    )

#    if model=='CanESM5' or model=='GISS-E2-1-G':
#        dirhist  = [x for x in dirhist if 'r1i1p1f1' in x]
    # experiment missing? skip model
    if len(path_hist_tas)==0:
        print(model + ' not provided historical tas')
        continue
    historical[model]['tas'] = np.zeros((165))
    nens = 0
    for ens in path_hist_tas:
        print(ens)
        tas = np.loadtxt(ens)
        if tas.size >= 165:
            historical[model]['tas']  = historical[model]['tas'] + tas[:165]
            nens = nens + 1
    if nens == 0:
        continue
    historical[model]['tas'] = historical[model]['tas'] / nens
    nyears[model]  = 165
    historical[model]['1951-1980'] = np.mean(historical[model]['tas'][101:131]) - np.mean(historical[model]['tas'][0:51])
    historical[model]['1961-1990'] = np.mean(historical[model]['tas'][111:141]) - np.mean(historical[model]['tas'][0:51])
    historical[model]['1995-2014'] = np.mean(historical[model]['tas'][145:165]) - np.mean(historical[model]['tas'][0:51])
    # if we get this far, things have worked out well
    accepted_models.append(model)

# %%
len(accepted_models)
#nyears

# %%
cw_temp = np.loadtxt('../data_input/CW.txt')
blratio = np.loadtxt('../data_input/cmip5_data_2019.txt')[5,:]
cowtan = cw_temp[:,1] - np.mean(cw_temp[:51,1])
blratio  = np.concatenate((np.ones(11), blratio))
Tobs = blratio * cowtan
#Tobs=cowtan
print(np.mean(Tobs[111:141]))
print(np.mean(Tobs[101:131]))

# %%
sixtyoneninety=np.ones(len(accepted_models))*np.nan
fiftyoneeighty=np.ones(len(accepted_models))*np.nan
ninetyfivefourteen = np.ones(len(accepted_models))*np.nan
full=np.ones((165, len(accepted_models)))
for i, model in enumerate(accepted_models):
    full[:,i] = historical[model]['tas'][:165] - np.mean(historical[model]['tas'][0:51])
    pl.plot(np.arange(1850, 1850+nyears[model]), historical[model]['tas'] - np.mean(historical[model]['tas'][0:51]))
    sixtyoneninety[i] = historical[model]['1961-1990']
    fiftyoneeighty[i] = historical[model]['1951-1980']
    ninetyfivefourteen[i] = historical[model]['1995-2014']
pl.plot(np.arange(1850, 2020), Tobs, color='k', lw=2)

# %%
fig, ax=pl.subplots()#figsize=(9.5/2.54,9.5/2.54))
ax.fill_between(np.arange(1850.5,2015), np.mean(full,axis=1)-np.std(full, axis=1), np.mean(full,axis=1)+np.std(full,axis=1), color='green', alpha=0.5)
ax.plot(np.arange(1850.5,2015), np.mean(full, axis=1), color='green', label='CMIP6 historical')
ax.fill_between(np.arange(1850.5,2015), Tobs[:-5]-cw_temp[:-5,2], Tobs[:-5]+cw_temp[:-5,2], color='k', alpha=0.5)
ax.plot(np.arange(1850.5,2015), Tobs[:-5], color='k', label='Reconstructed GSAT')
ax.set_xlim(1850,2015)
ax.set_ylim(-0.4, 1.35)
ax.legend(loc='upper left')
ax.set_ylabel('Temperature anomaly with respect to 1850-1900, $^{\circ}$C')
ax.set_title('CMIP6 simulated and observed warming')
pl.tight_layout()
pl.savefig('../figures/figureS7.png', dpi=300)
pl.savefig('../figures/figureS7.pdf')

# %%
print(np.mean(sixtyoneninety))
print(np.mean(fiftyoneeighty))

# %%
print(np.std(sixtyoneninety))
print(np.std(fiftyoneeighty))

# %%
# cowtan and way uncertainty from 1850-1900 to 1961-90 (one sigma)
np.sqrt(np.sqrt(np.sum(cw_temp[:51,2]**2)/51)**2 + np.sqrt(np.sum(cw_temp[111:141,2]**2)/30)**2)

# %%
for model in ['CanESM5','E3SM-1-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','IPSL-CM6A-LR',
             'MIROC6','MRI-ESM2-0','NorESM2-LM','UKESM1-0-LL']:
    pl.plot(historical[model]['tas'][95:121]-historical[model]['tas'][95])

# %%
for model in ['CanESM5','E3SM-1-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','IPSL-CM6A-LR',
             'MIROC6','MRI-ESM2-0','NorESM2-LM','UKESM1-0-LL']:
    print(model, historical[model]['1995-2014']-historical[model]['1951-1980'])

# %%
st.linregress(np.arange(11), Tobs[159:])

# %%
