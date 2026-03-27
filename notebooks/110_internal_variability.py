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
from statsmodels.tsa.stattools import acovf
import glob
import matplotlib.pyplot as pl
import warnings
from tqdm.auto import tqdm
from scipy.interpolate import interp1d
warnings.simplefilter('ignore')

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

# %%
models = [
'ACCESS-CM2',    'CAMS-CSM1-0','CESM2-WACCM-FV2','E3SM-1-1-ECA','GISS-E2-1-G-CC', 'INM-CM4-8',    'MPI-ESM1-2-LR','NorESM2-MM',
'ACCESS-ESM1-5', 'CanESM5',    'CNRM-CM6-1',     'FGOALS-f3-L', 'GISS-E2-1-H',    'INM-CM5-0',    'MRI-ESM2-0',   'SAM0-UNICON',
'AWI-CM-1-1-MR', 'CAS-ESM2-0', 'CNRM-CM6-1-HR',  'FGOALS-g3',   'GISS-E2-2-G',    'IPSL-CM6A-LR', 'NESM3',        'UKESM1-0-LL',
'AWI-ESM-1-1-LR','CESM2',      'CNRM-ESM2-1',    'GFDL-CM4',    'HadGEM3-GC31-LL','MIROC6',       'NorCPM1',
'BCC-CSM2-MR',   'CESM2-FV2',  'E3SM-1-0',       'GFDL-ESM4',   'HadGEM3-GC31-MM','MIROC-ES2L',   'NorESM1-F',
'BCC-ESM1',      'CESM2-WACCM','E3SM-1-1',       'GISS-E2-1-G', 'IITM-ESM',       'MPI-ESM1-2-HR','NorESM2-LM',
]

models = [
'ACCESS-CM2',     'CanESM5-CanOE',   'CNRM-CM6-1',    'EC-Earth3-LR',     'GISS-E2-1-G',     'INM-CM5-0',       'NESM3',
'ACCESS-ESM1-5',  'CAS-ESM2-0',      'CNRM-CM6-1-HR', 'EC-Earth3-Veg',    'GISS-E2-1-G-CC',  'IPSL-CM6A-LR',    'NorCPM1',
'AWI-CM-1-1-MR',  'CESM2',           'CNRM-ESM2-1',   'EC-Earth3-Veg-LR', 'GISS-E2-1-H',     'MIROC6',          'NorESM1-F',
'AWI-ESM-1-1-LR', 'CESM2-FV2',                        'FGOALS-f3-L',      'GISS-E2-2-G',     'MIROC-ES2L',      'NorESM2-LM',
'BCC-CSM2-MR',    'CESM2-WACCM',     'E3SM-1-0',      'FGOALS-g3',        'HadGEM3-GC31-LL', 'MPI-ESM-1-2-HAM', 'NorESM2-MM',
'BCC-ESM1',       'CESM2-WACCM-FV2', 'E3SM-1-1',      'FIO-ESM-2-0',      'HadGEM3-GC31-MM', 'MPI-ESM1-2-HR',   'SAM0-UNICON',
'CAMS-CSM1-0',    'CIESM',           'E3SM-1-1-ECA',  'GFDL-CM4',         'IITM-ESM',        'MPI-ESM1-2-LR',   'TaiESM1',
'CanESM5',        'CMCC-CM2-SR5',    'EC-Earth3',     'GFDL-ESM4',        'INM-CM4-8',       'MRI-ESM2-0',      'UKESM1-0-LL'
]

# %%
piControl = {}
accepted_models = []
nyears = {}

# Models only provided r1 in any seriousness. GISS and CanESM have a couple of variants, we'll stick to r1 for consistency.
# All of this data has been pre-processed into annual mean text files. It's way too big to put on GitHub.
for model in models:
    piControl[model] = {}
    dirpiC  = glob.glob('../data_input/cmip6/%s/piControl/r1i*/' % model)
    if model=='CanESM5' or model=='GISS-E2-1-G':
        dirpiC  = [x for x in dirpiC if 'r1i1p1f1' in x]
    # experiment missing? skip model
    if len(dirpiC)==0:
        continue
    dirpiC  = dirpiC[0]
    piControl[model]['tas']  = np.loadtxt(dirpiC + 'tas.txt')
    nyears[model]  = len(piControl[model]['tas'])
    slope, intercept, _, _, _ = st.linregress(np.arange(nyears[model]), piControl[model]['tas'])
    piControl[model]['tas_driftcorrected'] = piControl[model]['tas'] - intercept - np.arange(nyears[model]) * slope
    # if we get this far, things have worked out well
    accepted_models.append(model)

# %%
len(accepted_models)
#nyears

# %%
fig, ax = pl.subplots(10,5, figsize = (19/2.54,25/2.54))
for i, model in enumerate(sorted(accepted_models)):
    ax[i//5,i%5].plot(piControl[model]['tas_driftcorrected'], color='k')
    ax[i//5,i%5].set_title(model, fontsize=7)
    #ax[i//5,i%5].xaxis.set_label_coords(0.5, -0.12)
    #ax[i//5,i%5].set_ylabel('T (K)')
    #ax[i//5,i%5].yaxis.set_label_coords(-0.035, 0.9)
    ax[i//5,i%5].set_xlim(0,nyears[model])#(0,nyears[model])
    ax[i//5,i%5].set_ylim(-0.5,0.5)#(0,np.max(np.concatenate((tg[model],T[model])))))
    if i%5!=0:
        ax[i//5,i%5].set_yticklabels([])
ax[0,0].set_ylabel('T (K)')
ax[1,0].set_ylabel('T (K)')
ax[2,0].set_ylabel('T (K)')
ax[3,0].set_ylabel('T (K)')
ax[4,0].set_ylabel('T (K)')
ax[5,0].set_ylabel('T (K)')
ax[6,0].set_ylabel('T (K)')
ax[7,0].set_ylabel('T (K)')
ax[8,0].set_ylabel('T (K)')
ax[9,0].set_ylabel('T (K)')
#ax[8,0].set_xlabel('year')
#ax[8,1].set_xlabel('year')
#ax[8,2].set_xlabel('year')
#ax[8,3].set_xlabel('year')
#ax[7,4].set_xlabel('year')
ax[9,4].axis('off')
fig.tight_layout();
pl.savefig('../figures/figureS5.png', dpi=300)
pl.savefig('../figures/figureS5.pdf')


# %%
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]

fig, ax = pl.subplots()
for i, model in enumerate(accepted_models):
    y = autocorr(piControl[model]['tas_driftcorrected'])
    ax.plot(y / float(y.max()))
ax.set_xlim(0,276)

# %%
for i, model in enumerate(accepted_models):
    print(model, np.std(piControl[model]['tas_driftcorrected']))

# %%
acov  = {}
for model in accepted_models:
    acov[model] = np.zeros(276)
    nyears = len(piControl[model]['tas_driftcorrected'])
    if nyears<276:
        acov[model][:nyears] = acovf(piControl[model]['tas_driftcorrected'])
    else:
        acov[model] = acovf(piControl[model]['tas_driftcorrected'])[:276]

# %%
acm = {}
for model in accepted_models:
    for i in range(1, 276):
        acm[model] = np.zeros((276, 276))
        acm[model] = acm[model] + np.diag(acov[model][i]*np.ones(276-i), i) + np.diag(acov[model][i]*np.ones(276-i), -i)
    acm[model] = acm[model] + np.diag(acov[model][0]*np.ones(276))

# %%
acov = {}
acm = {}
nyears = {}
for model in tqdm(accepted_models):
    ac = acovf(piControl[model]['tas_driftcorrected'])
    nyears[model] = len(piControl[model]['tas_driftcorrected'])
    if nyears[model]<276:
        acov[model] = np.zeros(276)
        acov[model][:nyears[model]] = ac
        nyears[model] = 276
    else:
        acov[model] = ac
    acm[model] = np.zeros((nyears[model], nyears[model]))
    for i in range(1, nyears[model]):
        acm[model] = acm[model] + np.diag(acov[model][i]*np.ones(nyears[model]-i), i) + np.diag(acov[model][i]*np.ones(nyears[model]-i), -i)
    acm[model] = acm[model] + np.diag(acov[model][0]*np.ones(nyears[model]))

# %%
#for model in models:
x = st.multivariate_normal.rvs(cov=acm['CNRM-ESM2-1'], random_state=10)
y = st.multivariate_normal.rvs(cov=acm['CNRM-ESM2-1'], random_state=11)
z = st.multivariate_normal.rvs(cov=acm['CNRM-ESM2-1'], random_state=12)

pl.plot(x[:276])
pl.plot(y[:276])
pl.plot(z[:276])


# %%
np.random.seed(seed=360185)
samples = 100000
model_choices = np.random.randint(0, high=len(accepted_models), size=samples)
intvar = np.zeros((276, samples))
for i in tqdm(range(samples)):
    intvar[:, i] = st.multivariate_normal.rvs(cov=acm[accepted_models[model_choices[i]]], random_state=98426+i*9)[:276]

# %%
os.makedirs('../data_output/piControl', exist_ok=True)

# %%
np.savetxt('../data_output/piControl/internal_variability_piControl.txt', intvar)

# %%
