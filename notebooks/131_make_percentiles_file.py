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
import h5py
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from netCDF4 import Dataset
from tqdm.auto import tqdm
import wquantiles

# %%
expts = ['CMIP6-constrained','CanESM5','E3SM','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','HadGEM3-GC31-LL','IPSL-CM6A-LR','MIROC6','MRI-ESM2-0','NorESM2-LM','Oslo-CTM3','UKESM1-0-LL']
#expts_all = ['CMIP6-SSP1-1.9','CMIP6-SSP2-4.5','CMIP6-SSP3-7.0','CanESM5','E3SM','GFDL-CM4','GISS-E2-1-G','HadGEM3-GC31-LL','MIROC6','NorESM2-LM','Lund','AR5']

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
ERFari = load_dict_from_hdf5('../data_output/results/ERFari.h5')
ERFaci = load_dict_from_hdf5('../data_output/results/ERFaci.h5')
temp   = load_dict_from_hdf5('../data_output/results/temp.h5')
ks     = load_dict_from_hdf5('../data_output/results/knutti_score.h5')
ohc    = load_dict_from_hdf5('../data_output/results/ohc.h5')
hflux  = load_dict_from_hdf5('../data_output/results/hflux.h5')

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
colors = {
    'ECLIPSE-constrained'  : '0.6',
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

# %%
nc = Dataset('../data_input/GCOS_all_heat_content_1960-2018_ZJ_v22062020.nc')
#print(nc.variables)
ohctop = nc.variables['ohc_0-2000m'][:]
ohcbot = nc.variables['ohc_below_2000m'][:]
atmosh = nc.variables['atmospheric_heat_content'][:]
cryosh = nc.variables['energy_cryosphere'][:]
landhc = nc.variables['ground_heat_content'][:]
ohctopu = nc.variables['ohc_0-2000m_uncertainty'][:]
ohcbotu = nc.variables['ohc_below_2000m_uncertainty'][:]
atmoshu = nc.variables['atmospheric_heat_content_uncertainty'][:]
cryoshu = nc.variables['energy_cryosphere_uncertainty'][:]
landhcu = nc.variables['ground_heat_content_uncertainty'][:]
nc.close()
cryosh[-1] = cryosh[-2]  # nan for 2017-2018, assume no change
cryoshu[-1] = cryoshu[-2]
OHCobs = (ohctop+ohcbot+atmosh+cryosh+landhc)-(ohctop+ohcbot+atmosh+cryosh+landhc)[11]

# %%
fig, ax = pl.subplots()
for expt in expts:
    ax.plot(np.arange(1750.5,2020), 10*np.nansum((ohc[expt]-ohc[expt][221,:])*ks['multi'][expt], axis=1), color=colors[expt])
ax.set_title('Best estimate aerosol forcing')
ax.set_ylabel('Ocean heat uptake relative to 1960, ZJ')
ax.plot(np.arange(1960.5,2019), (OHCobs-OHCobs[11]), color=colors['observations'])
ax.set_xlim(1960,2019)
ax.set_ylim(-50,400)

# %%
geoff_sample_df = pd.read_csv('../data_output/geoff_sample.csv')
geoff_sample_df

# %%
# Throw TCR into the mix
samples = 100000
tcr = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'] + geoff_sample_df[:samples]['eff']*geoff_sample_df[:samples]['gamma_2l'])
ecs = geoff_sample_df['q4x'][:samples]/2/(geoff_sample_df[:samples]['lamg'])


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
                pc[expt][constraint][metric][perc] = np.zeros(270)
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
        for year in range(270):
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
save_dict_to_hdf5(pc, '../data_output/results/pc.h5')

# %%
