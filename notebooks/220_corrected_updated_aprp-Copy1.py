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
# # Use updated APRP decomposition and 13 models
#
# The original APRP code I did was wrong (Zelinka et al. 2023), so use the updated fits that I did in https://github.com/chrisroadmap/cmip6-aerosol-forcing.

# %%
import pandas as pd
import glob
import os
import matplotlib.pyplot as pl
from scipy.signal import savgol_filter
import numpy as np
from scipy.optimize import curve_fit
import json
import scipy

# %%
models = [
    'CanESM5', 
    'CNRM-CM6-1',
    'E3SM-2-0', 
    'GFDL-CM4', 
    'GFDL-ESM4', 
    'GISS-E2-1-G',
    'HadGEM3-GC31-LL',
    'IPSL-CM6A-LR',
    'MIROC6',
    'MPI-ESM-1-2-HAM',
    'MRI-ESM2-0',
    'NorESM2-LM',
    'UKESM1-0-LL'
]

colors = {
    "CanESM5": "red",
    "CNRM-CM6-1": "orangered",
    "E3SM-2-0": "darkorange",
    "GFDL-ESM4": "yellowgreen",
    "GFDL-CM4": "yellow",
    "GISS-E2-1-G": "green",
    "HadGEM3-GC31-LL": "turquoise",
    "IPSL-CM6A-LR": "teal",
    "MIROC6": "blue",
    "MPI-ESM-1-2-HAM": "darkslateblue",
    "MRI-ESM2-0": "blueviolet",
    "NorESM2-LM": "purple",
    "UKESM1-0-LL": "crimson",
}

# %%
models_runs = {}

for path in glob.glob(os.path.join(*os.path.normpath('../data_input/corrected_aprp/*/*/mean/*.csv').split(os.sep))):
    model = path.split(os.sep)[3]
    if model not in models_runs:
        models_runs[model] = []
    run = path.split(os.sep)[4]
    models_runs[model].append(run)

ari = {}
aci = {}
erf = {}
years = {}

for model in models_runs:
    nruns = 0
    for run in models_runs[model]:
        path = os.path.join('..', 'data_input', 'corrected_aprp', model, run, 'mean', f'{model}_{run}_aerosol_forcing.csv')
        df = pd.read_csv(path, index_col=0)
        if nruns == 0:
            erf_temp = df['ERF'].values.squeeze()
            ari_temp = df["ERFari"].values.squeeze()
            aci_temp = df["ERFaci"].values.squeeze()
        else:
            erf_temp = erf_temp + df['ERF'].values.squeeze()
            ari_temp = ari_temp + df["ERFari"].values.squeeze()
            aci_temp = aci_temp + df["ERFaci"].values.squeeze()
        years[model] = df.index
        nruns = nruns + 1
    erf[model] = erf_temp / nruns
    ari[model] = ari_temp / nruns
    aci[model] = aci_temp / nruns
    #print(os.path.normpath('../data_input/corrected_aprp/*/*/mean/*.csv').split(os.sep))
    #files.append(os.path.split(path)[-1])

# %%
fig, ax = pl.subplots(figsize=(19/2.54,9/2.54))
for model in models:
    ax.plot(years[model], erf[model], '.', lw=0, alpha=0.3, color=colors[model])
    ax.plot(years[model], savgol_filter(np.concatenate((np.zeros(10), erf[model], erf[model][-1]*np.ones(10))), 11, 1)[:-20], color=colors[model], label=model)
ax.set_xlim(1850,2025)
ax.set_ylim(-2.2,0.7)
ax.plot(2014.5,-1.01,'ko')
ax.plot((2014.5,2014.5),(-0.63,-1.37),'k',lw=3,solid_joinstyle='bevel')
ax.legend(ncol=2, frameon=False)
ax.set_title('CMIP6 aerosol effective radiative forcing with respect to 1850: historical + SSP2-4.5')
ax.set_ylabel('W m$^{-2}$')
fig.tight_layout()
# pl.savefig('../figures/figure1.png', dpi=300)
# pl.savefig('../figures/figure1.pdf')

# %%
emis_df = pd.read_csv('../data_input/rcmip/rcmip-emissions-annual-means-v5-1-0.csv')

bc = (
    emis_df.loc[
        (emis_df["Scenario"] == "ssp245")
        & (emis_df["Region"] == "World")
        & (emis_df["Variable"] == "Emissions|BC"),
        "1750":"2100",
    ]
    .interpolate(axis=1)
    .squeeze()
    .values
)

oc = (
    emis_df.loc[
        (emis_df["Scenario"] == "ssp245")
        & (emis_df["Region"] == "World")
        & (emis_df["Variable"] == "Emissions|OC"),
        "1750":"2100",
    ]
    .interpolate(axis=1)
    .squeeze()
    .values
)

so2 = (
    emis_df.loc[
        (emis_df["Scenario"] == "ssp245")
        & (emis_df["Region"] == "World")
        & (emis_df["Variable"] == "Emissions|Sulfur"),
        "1750":"2100",
    ]
    .interpolate(axis=1)
    .squeeze()
    .values
)


# %%
def ari_linear(x, a0, a1, a2):
    ari = x[0] * a0 + x[1] * a1 + x[2] * a2
    ari_1850 = so2[100] * a0 + bc[100] * a1 + oc[100] * a2
    return ari - ari_1850

def aci_log(x, beta, n0, n1, n2):
    aci = beta * np.log(1 + x[0] * n0 + x[1] * n1 + x[2] * n2)
    aci_1850 = beta * np.log(1 + so2[100] * n0 + bc[100] * n1 + oc[100] * n2)
    return aci - aci_1850


# %%
ari_fits = {}
aci_fits = {}

for model in models:
    print(model)
    
    ist = int(np.floor(years[model][0] - 1750))
    ien = int(np.ceil(years[model][-1] - 1750 + 1))
    ari_fits[model], cov = curve_fit(
        ari_linear,
        [so2[ist:ien], bc[ist:ien], oc[ist:ien]],
        ari[model],
    )
    
    ist = int(np.floor(years[model][0] - 1750))
    ien = int(np.ceil(years[model][-1] - 1750 + 1))
    aci_fits[model], cov = curve_fit(
        aci_log,
        [so2[ist:ien], bc[ist:ien], oc[ist:ien]],
        aci[model],
        bounds=((-np.inf, 0, 0, 0), (0, np.inf, np.inf, np.inf)),
        max_nfev=10000,
    )

# %%
fig, ax = pl.subplots(4, 4, figsize=(24 / 2.54, 16 / 2.54), squeeze=False)
for imodel, model in enumerate(sorted(models, key=str.lower)):
    i = imodel // 4
    j = imodel % 4
    ax[i, j].plot(years[model], aci[model], color="k", ls="-", alpha=0.5, lw=1)
    ax[i, j].plot(
        np.arange(1750.5, 2101),
        aci_log([so2, bc, oc], *aci_fits[model]),
        color=colors[model],
        zorder=7,
        lw=1,
    )

    ax[i, j].set_xlim(1750, 2100)
    ax[i, j].set_ylim(-1.7, 0.5)
    ax[i, j].axhline(0, lw=0.5, ls=":", color="k")
    ax[i, j].fill_between(
        np.arange(1850, 2015), -10, 10, color="#e0e0e0", zorder=-20
    )
    ax[i, j].get_xticklabels()[-1].set_ha("right")
    if model == "HadGEM3-GC31-LL":
        modlab = "HadGEM3"
    elif model == "MPI-ESM-1-2-HAM":
        modlab = "MPI-ESM1-2"
    else:
        modlab = model
    ax[i, j].text(
        0.03, 0.05, modlab, transform=ax[i, j].transAxes, fontweight="bold"
    )

ax[0, 0].set_ylabel("W m$^{-2}$")
ax[1, 0].set_ylabel("W m$^{-2}$")
ax[2, 0].set_ylabel("W m$^{-2}$")
ax[3, 0].set_ylabel("W m$^{-2}$")
ax[3, 1].axis("off")
ax[3, 2].axis("off")
ax[3, 3].axis("off")

pl.suptitle("Aerosol-cloud interactions parameterisations")

fig.tight_layout()

# %%
fig, ax = pl.subplots(4, 4, figsize=(24 / 2.54, 16 / 2.54), squeeze=False)
for imodel, model in enumerate(sorted(models, key=str.lower)):
    i = imodel // 4
    j = imodel % 4
    ax[i, j].plot(years[model], ari[model], color="k", ls="-", alpha=0.5, lw=1)
    ax[i, j].plot(
        np.arange(1750.5, 2101),
        ari_linear([so2, bc, oc], *ari_fits[model]),
        color=colors[model],
        zorder=7,
        lw=1,
    )

    ax[i, j].set_xlim(1750, 2100)
    ax[i, j].set_ylim(-1, 0.5)
    ax[i, j].axhline(0, lw=0.5, ls=":", color="k")
    ax[i, j].fill_between(
        np.arange(1850, 2015), -10, 10, color="#e0e0e0", zorder=-20
    )
    ax[i, j].get_xticklabels()[-1].set_ha("right")
    if model == "HadGEM3-GC31-LL":
        modlab = "HadGEM3"
    elif model == "MPI-ESM-1-2-HAM":
        modlab = "MPI-ESM1-2"
    else:
        modlab = model
    ax[i, j].text(
        0.03, 0.05, modlab, transform=ax[i, j].transAxes, fontweight="bold"
    )

ax[0, 0].set_ylabel("W m$^{-2}$")
ax[1, 0].set_ylabel("W m$^{-2}$")
ax[2, 0].set_ylabel("W m$^{-2}$")
ax[3, 0].set_ylabel("W m$^{-2}$")
ax[3, 1].axis("off")
ax[3, 2].axis("off")
ax[3, 3].axis("off")

pl.suptitle("Aerosol-radiation interactions parameterisations")

fig.tight_layout()

# %%
aci_fits

# %%
coeff = {}
for model in aci_fits:
    coeff[model] = {}
    coeff[model]['ERFaci'] = {
        'beta': aci_fits[model][0],
        'n0': aci_fits[model][1],
        'n1': aci_fits[model][2],
        'n2': aci_fits[model][3]
    }
    coeff[model]['ERFari'] = {
        'SO2': aci_fits[model][0],
        'BC': aci_fits[model][1],
        'OC': aci_fits[model][2]
    }

# %%
so2_samp_ari = np.zeros(len(models))
bc_samp_ari = np.zeros(len(models))
oc_samp_ari = np.zeros(len(models))
log_n0_samp = np.zeros(len(models))
log_n1_samp = np.zeros(len(models))
log_n2_samp = np.zeros(len(models))
beta_samp = np.zeros(len(models))

for im, model in enumerate(models):
    so2_samp_ari[im] = coeff[model]['ERFari']['SO2']
    bc_samp_ari[im] = coeff[model]['ERFari']['BC']
    oc_samp_ari[im] = coeff[model]['ERFari']['OC']
    beta_samp[im] = (coeff[model]['ERFaci']['beta'])
    log_n0_samp[im] = np.log(coeff[model]['ERFaci']['n0'])
    log_n1_samp[im] = np.log(coeff[model]['ERFaci']['n1'])
    log_n2_samp[im] = np.log(coeff[model]['ERFaci']['n2'])

# %%
with open("../data_output/corrected_cmip6_aerosol_coefficients.json", "w") as write_file:
    json.dump(coeff, write_file, indent=4)

# %%
samples = 100000
kde = scipy.stats.gaussian_kde(
    [log_n0_samp, log_n1_samp, log_n2_samp], bw_method=0.1
)
aci_sample = kde.resample(size=samples * 1, seed=63648708)

# %%
aci_sample

# %%
kde = scipy.stats.gaussian_kde([so2_samp_ari, bc_samp_ari, oc_samp_ari])
ari_sample=kde.resample(size=samples, seed=685534562)

# %%
ari_sample

# %%
df = pd.DataFrame(aci_sample.T, columns=['log(n_SO2)', 'log(n_BC)', 'log(n_OC)'])
df.to_csv('../data_output/corrected_ERFaci_samples.csv', index=False)

df = pd.DataFrame(ari_sample.T, columns=['SO2', 'BC', 'OC'])
df.to_csv('../data_output/corrected_ERFari_samples.csv', index=False)

# %%
