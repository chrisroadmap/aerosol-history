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
import matplotlib.gridspec as gridspec
import matplotlib.font_manager
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D
import matplotlib.pyplot as pl

import numpy as np
import scipy.stats as st
import pandas as pd
import os
import wquantiles
from tqdm.auto import tqdm
from netCDF4 import Dataset
import h5py


# %%
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
            ans[key] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans


# %%
# load large datafiles calculated previously
ERFari = load_dict_from_hdf5('../data_output/results/igcc_2022/ERFari.h5')
ERFaci = load_dict_from_hdf5('../data_output/results/igcc_2022/ERFaci.h5')
temp   = load_dict_from_hdf5('../data_output/results/igcc_2022/temp.h5')
ks     = load_dict_from_hdf5('../data_output/results/igcc_2022/knutti_score.h5')
ohc    = load_dict_from_hdf5('../data_output/results/igcc_2022/ohc.h5')
hflux  = load_dict_from_hdf5('../data_output/results/igcc_2022/hflux.h5')

# %%
samples = 100000
geoff_sample_df = pd.read_csv('../data_output/geoff_sample.csv')
geoff_sample_df

# %%
# Throw TCR into the mix
samples = 100000
tcr = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'] + geoff_sample_df[:samples]['eff']*geoff_sample_df[:samples]['gamma_2l'])
ecs = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'])

# %%
expts = ['CMIP6-constrained']


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
intvar = np.loadtxt('../data_output/piControl/internal_variability_piControl.txt')

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
for constraint in ['temp','ohc','multi']:
    print(constraint)
    print('---------')
    for expt in expts:
        print(expt, constraint, 'ECS', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ECS']['5'], pc[expt][constraint]['ECS']['16'], pc[expt][constraint]['ECS']['50'], np.sum(ecs*ks[constraint][expt]), pc[expt][constraint]['ECS']['84'], pc[expt][constraint]['ECS']['95']))
        print(expt, constraint, 'TCR', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['TCR']['5'], pc[expt][constraint]['TCR']['16'], pc[expt][constraint]['TCR']['50'], np.sum(tcr*ks[constraint][expt]), pc[expt][constraint]['TCR']['84'], pc[expt][constraint]['TCR']['95']))
        print(expt, constraint, 'ERFaer2022', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFaer']['5'][272], pc[expt][constraint]['ERFaer']['16'][272], pc[expt][constraint]['ERFaer']['50'][272], np.sum((ERFari[expt][272] + ERFaci[expt][272])*ks[constraint][expt]), pc[expt][constraint]['ERFaer']['84'][272], pc[expt][constraint]['ERFaer']['95'][272]))
        print(expt, constraint, 'ERFari2022', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFari']['5'][272], pc[expt][constraint]['ERFari']['16'][272], pc[expt][constraint]['ERFari']['50'][272], np.sum(ERFari[expt][272]*ks[constraint][expt]), pc[expt][constraint]['ERFari']['84'][272], pc[expt][constraint]['ERFari']['95'][272]))
        print(expt, constraint, 'ERFaci2022', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFaci']['5'][272], pc[expt][constraint]['ERFaci']['16'][272], pc[expt][constraint]['ERFaci']['50'][272], np.sum(ERFaci[expt][272]*ks[constraint][expt]), pc[expt][constraint]['ERFaci']['84'][272], pc[expt][constraint]['ERFaci']['95'][272]))
        print(expt, constraint, 'ERFaer2019', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFaer']['5'][269], pc[expt][constraint]['ERFaer']['16'][269], pc[expt][constraint]['ERFaer']['50'][269], np.sum((ERFari[expt][269] + ERFaci[expt][269])*ks[constraint][expt]), pc[expt][constraint]['ERFaer']['84'][269], pc[expt][constraint]['ERFaer']['95'][269]))
        print(expt, constraint, 'ERFari2019', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFari']['5'][269], pc[expt][constraint]['ERFari']['16'][269], pc[expt][constraint]['ERFari']['50'][269], np.sum(ERFari[expt][269]*ks[constraint][expt]), pc[expt][constraint]['ERFari']['84'][269], pc[expt][constraint]['ERFari']['95'][269]))
        print(expt, constraint, 'ERFaci2019', '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f' % (pc[expt][constraint]['ERFaci']['5'][269], pc[expt][constraint]['ERFaci']['16'][269], pc[expt][constraint]['ERFaci']['50'][269], np.sum(ERFaci[expt][269]*ks[constraint][expt]), pc[expt][constraint]['ERFaci']['84'][269], pc[expt][constraint]['ERFaci']['95'][269]))

# %%
colors = {
    'CMIP6-constrained' : '0.3',
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
    'observations'   : 'black',
    'Oslo-CTM3'      : 'pink',
}

ls = {
    'CMIP6-constrained'  : '-',
    'CanESM5'        : '-',
    'E3SM'           : '-',
    'GFDL-ESM4'      : '-',
    'GFDL-CM4'       : '-',
    'GISS-E2-1-G'    : '-',
    'HadGEM3-GC31-LL': '-',
    'IPSL-CM6A-LR'   : '-',
    'MIROC6'         : '-',
    'MRI-ESM2-0'     : '-',
    'NorESM2-LM'     : '-',
    'UKESM1-0-LL'    : '-',
    'Oslo-CTM3'      : '-',
    'observations'   : '-',
}

# %%
df_eei = pd.read_csv('../data_input/IGCC2022_earth_energy_imbalance.csv', index_col=0)
OHCobs = (df_eei.loc[:2020.5, 'Total'] - df_eei.loc[1971.5, 'Total']).values

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
fig, ax = pl.subplots()
for expt in expts:
    ax.plot(np.arange(1750.5,2023), 10*np.nansum((ohc[expt]-ohc[expt][221,:])*ks['multi'][expt], axis=1), color=colors[expt], label=expt)
ax.set_title('Best estimate aerosol forcing')
ax.set_ylabel('Ocean heat uptake relative to 1971, ZJ')
ax.plot(np.arange(1971.5,2021), (OHCobs-OHCobs[0]), color=colors['observations'], label='obs., IGCC 2022')
ax.set_xlim(1971,2021)
ax.set_ylim(-50,500)
ax.legend()

# %%
# Temperature (GMST) observations: IGCC 2022
tempobs = pd.read_csv('../data_input/IGCC2022_annual_averages.csv')['gmst'].values
blratio = 1
years  = np.arange(1850.5, 2023)
Tobs = blratio * tempobs

# %%
#fig, ax = pl.subplots(3,2,figsize=(16/2.54,23/2.54))
fig = pl.figure(figsize=(19/2.54,14/2.54))
gs = gridspec.GridSpec(ncols=6, nrows=2)
axtmp = fig.add_subplot(gs[0,4:6])
axohc = fig.add_subplot(gs[1,4:6])
axari = fig.add_subplot(gs[0,0:2])
axaci = fig.add_subplot(gs[0,2:4])
axaer = fig.add_subplot(gs[1,0:4])

expt='CMIP6-constrained'
axtmp.fill_between(np.arange(1750.5,2023), pc[expt]['multi']['GSAT']['5'], pc[expt]['multi']['GSAT']['95'], color=colors[expt], alpha=0.3, lw=0)
axaer.fill_between(np.arange(1750.5,2023), pc[expt]['multi']['ERFaer']['5'], pc[expt]['multi']['ERFaer']['95'], color=colors[expt], alpha=0.3, lw=0, label='CMIP6-constrained 5-95% range')
axari.fill_between(np.arange(1750.5,2023), pc[expt]['multi']['ERFari']['5'], pc[expt]['multi']['ERFari']['95'], color=colors[expt], alpha=0.3, lw=0)
axaci.fill_between(np.arange(1750.5,2023), pc[expt]['multi']['ERFaci']['5'], pc[expt]['multi']['ERFaci']['95'], color=colors[expt], alpha=0.3, lw=0)
axohc.fill_between(np.arange(1750.5,2023), 10*(pc[expt]['multi']['OHC']['5']-pc[expt]['multi']['OHC']['5'][221]), 10*(pc[expt]['multi']['OHC']['95']-pc[expt]['multi']['OHC']['95'][221]), color=colors[expt], alpha=0.3, lw=0)

for expt in expts:
    if expt=='CMIP6-constrained':
        lw=2
        zorder=10
        axtmp.plot(np.arange(1750.5,2023), np.nansum((temp[expt] + intvar[:273,:samples])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], color=colors[expt], zorder=zorder)
        axaer.plot(np.arange(1750.5,2023), np.nansum((ERFari[expt][:273]+ERFaci[expt][:273])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axari.plot(np.arange(1750.5,2023), np.nansum((ERFari[expt][:273])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axaci.plot(np.arange(1750.5,2023), np.nansum((ERFaci[expt][:273])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axohc.plot(np.arange(1750.5,2023), 10*np.nansum((ohc[expt]-ohc[expt][221,:])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], color=colors[expt], zorder=zorder)

    else:
        lw=1
        zorder=1
        axtmp.plot(np.arange(1750.5,2023), np.nansum((temp[expt] + intvar[:273,:samples])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], color=colors[expt], zorder=zorder)
        axaer.plot(np.arange(1750.5,2023), savgol_filter(np.nansum((ERFari[expt][:273]+ERFaci[expt][:273])*ks['multi'][expt], axis=1), 11, 1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axari.plot(np.arange(1750.5,2023), savgol_filter(np.nansum((ERFari[expt][:273])*ks['multi'][expt], axis=1), 11, 1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axaci.plot(np.arange(1750.5,2023), savgol_filter(np.nansum((ERFaci[expt][:273])*ks['multi'][expt], axis=1), 11, 1), lw=lw, ls=ls[expt], label=expt, color=colors[expt], zorder=zorder)
        axohc.plot(np.arange(1750.5,2023), 10*np.nansum((ohc[expt]-ohc[expt][221,:])*ks['multi'][expt], axis=1), lw=lw, ls=ls[expt], color=colors[expt], zorder=zorder)

axtmp.plot(years, Tobs, color=colors['observations'], label='Reconstructed GSAT', lw=1)
axohc.plot(np.arange(1971.5,2021), (OHCobs-OHCobs[0]), color=colors['observations'], lw=1, label='Total Earth energy uptake\n(GCOS)')
axtmp.legend(frameon=False)
axohc.legend(frameon=False, loc='upper left')
axtmp.set_xlim(1750,2023)
axtmp.set_ylim(-0.5,1.7)
axtmp.set_ylabel('anomaly since 1850-1900 (K)');
axtmp.set_title('(d) Temperature')
axtmp.axhline(0, ls=':', color='k')
axohc.axhline(0, ls=':', color='k')
axohc.set_ylabel('anomaly since 1971 (ZJ)')
axohc.set_title('(e) Earth energy uptake')
axohc.set_ylim(-50,500)
axohc.set_xlim(1960,2023)
axaer.set_xlim(1750,2023)
axaer.set_ylim(-2.0,0.1)
axaer.set_ylabel('W m$^{-2}$');
axaer.set_title('(c) Total aerosol ERF')
axaer.axhline(0, ls=':', color='k')
axaer.legend(fontsize=8, frameon=False, ncol=2);
axari.set_xlim(1750,2023)
axari.set_ylim(-1.4,0.15)
axari.set_ylabel('W m$^{-2}$');
axari.set_title('(a) ERFari')
axari.axhline(0, ls=':', color='k')
axaci.set_xlim(1750,2023)
axaci.set_ylim(-1.4,0.15)
axaci.set_ylabel('W m$^{-2}$');
axaci.set_title('(b) ERFaci')
axaci.axhline(0, ls=':', color='k')
fig.tight_layout()
pl.savefig('../figures/figure6_igcc2022.png', dpi=300)
pl.savefig('../figures/figure6_igcc2022.pdf')

# %%
