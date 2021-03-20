import numpy as np
from netCDF4 import Dataset
import os, errno
import sys
import glob
import iris
import iris.coord_categorisation
import iris.analysis.cartography
from iris.analysis.maths import multiply, log, exp, apply_ufunc
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import sys
import warnings

from climateforcing.aprp import create_input, calc_aprp


warnings.simplefilter('ignore')


version = iris.__version__
if int(version.split('.')[0]) < 2:
    iris.FUTURE.clip_latitudes=True


def mkdir_p(path):
    """Verbatim from 
    https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    and unbelievably useful."""

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def make_deltas_e3sm(baserun, pertrun, model, expt, make_aprp):
    basedir = '/nfs/a65/pmcjs/E3SM/%s/Amon/' % baserun
    pertdir = '/nfs/a65/pmcjs/E3SM/%s/Amon/' % pertrun
    deltadir = '/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/deltas/'+pertrun+'/piClim-histaer/'
    mkdir_p(deltadir)
    
    # make deltas and ERF
    rsdt_base = iris.load(basedir+'rsdt/gr/v20200528/rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut/gr/v20200528/rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut/gr/v20200528/rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    if make_aprp:
        rlutcs_base = iris.load(basedir+'rlutcs/gr/v20200528/rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_base)
        rlutcs_base = rlutcs_base.concatenate()[0]
        rsds_base = iris.load(basedir+'rsds/gr/v20200528/rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_base)
        rsds_base = rsds_base.concatenate()[0]
        rsus_base = iris.load(basedir+'rsus/gr/v20200528/rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_base)
        rsus_base = rsus_base.concatenate()[0]
        rsutcs_base = iris.load(basedir+'rsutcs/gr/v20200528/rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_base)
        rsutcs_base = rsutcs_base.concatenate()[0]
        clt_base = iris.load(basedir+'clt/gr/v20200528/clt*.nc')
        equalise_attributes(clt_base)
        clt_base = clt_base.concatenate()[0]

    rsdt_pert = iris.load(pertdir+'rsdt/gr/v20200528/rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut/gr/v20200528/rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut/gr/v20200528/rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    if make_aprp:
        rsds_pert = iris.load(pertdir+'rsds/gr/v20200528/rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_pert)
        unify_time_units(rsds_pert)
        rsds_pert = rsds_pert.concatenate()[0]
        rsus_pert = iris.load(pertdir+'rsus/gr/v20200528/rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_pert)
        unify_time_units(rsus_pert)
        rsus_pert = rsus_pert.concatenate()[0]
        rlutcs_pert = iris.load(pertdir+'rlutcs/gr/v20200528/rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_pert)
        unify_time_units(rlutcs_pert)
        rlutcs_pert = rlutcs_pert.concatenate()[0]
        rsutcs_pert = iris.load(pertdir+'rsutcs/gr/v20200528/rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_pert)
        unify_time_units(rsutcs_pert)
        rsutcs_pert = rsutcs_pert.concatenate()[0]
        clt_pert = iris.load(pertdir+'clt/gr/v20200528/clt*.nc')
        equalise_attributes(clt_pert)
        unify_time_units(clt_pert)
        clt_pert = clt_pert.concatenate()[0]

    nmonths = rsdt_pert.data.shape[0]
    nmonthsbase = rsdt_base.data.shape[0]
    nyears = nmonths//12
    nyearsbase = nmonthsbase//12
    print(nmonths)

    erf_sw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_sw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))

    erf_sw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_sw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))

    erf_sw_data = (rsdt_pert.data-rsut_pert.data)-(rsdt_base.data-rsut_base.data)
    erf_lw_data = (-rlut_pert.data)-(-rlut_base.data)
    if make_aprp:
        erf_sw_cs_data = (rsdt_pert.data-rsutcs_pert.data)-(rsdt_base.data-rsutcs_base.data)
        erf_lw_cs_data = (-rlutcs_pert.data)-(-rlutcs_base.data)

    erf_sw = iris.cube.Cube(
        erf_sw_data,
        var_name = 'erf_sw',
        long_name = "Shortwave effective radiative forcing at TOA",
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_sw.coord('time').attributes = {}
    erf_sw.coord('latitude').attributes = {}
    erf_sw.coord('longitude').attributes = {}
    iris.save(erf_sw, deltadir+'erf_sw.nc')
    
    erf_lw = iris.cube.Cube(
        erf_lw_data,
        var_name = 'erf_lw',
        long_name = 'Longwave effective radiative forcing at TOA',
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_lw.coord('time').attributes = {}
    erf_lw.coord('latitude').attributes = {}
    erf_lw.coord('longitude').attributes = {}
    iris.save(erf_lw, deltadir+'erf_lw.nc')

    erf = erf_lw + erf_sw.data
    erf.var_name = 'erf'
    erf.long_name = 'Effective radiative forcing at TOA'
    iris.save(erf, deltadir+'erf.nc')

    if make_aprp:    
        erf_sw_cs = iris.cube.Cube(
            erf_sw_cs_data,
            var_name = 'erf_sw_cs',
            long_name = "Shortwave clear-sky effective radiative forcing at TOA",
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_sw_cs.coord('time').attributes = {}
        erf_sw_cs.coord('latitude').attributes = {}
        erf_sw_cs.coord('longitude').attributes = {}
        iris.save(erf_sw_cs, deltadir+'erf_sw_cs.nc')
        
        erf_lw_cs = iris.cube.Cube(
            erf_lw_cs_data,
            var_name = 'erf_lw_cs',
            long_name = 'Longwave clear-sky effective radiative forcing at TOA',
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_lw_cs.coord('time').attributes = {}
        erf_lw_cs.coord('latitude').attributes = {}
        erf_lw_cs.coord('longitude').attributes = {}
        iris.save(erf_lw_cs, deltadir+'erf_lw_cs.nc')
        
        erf_cs = erf_lw_cs + erf_sw_cs.data
        erf_cs.var_name = 'erf_cs'
        erf.long_name = 'Clear-sky effective radiative forcing at TOA'
        iris.save(erf_cs, deltadir+'erf_cs.nc')
    
    # easiest to do by working directly with the data
    adjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
    adjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'

    mkdir_p(adjalldir)
    mkdir_p(adjglobalmeanyearmeandir)

    # global means, statistics, etc
    variables = ['erf_sw', 'erf_lw', 'erf']
    if make_aprp:
        variables.extend(['erf_sw_cs', 'erf_lw_cs', 'erf_cs'])
    for var in variables:
        cube = iris.load_cube(deltadir+var+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, adjglobalmeanyearmeandir+var+'.nc')
    return


def make_deltas_aerchemmip(baserun, pertrun, model, expt, make_aprp, make_od550aer, has_rfmip):
    basedir = '/nfs/a65/pmcjs/CMIP6/AerChemMIP/'+model+'/'+baserun+'/histSST-piAer/'
    pertdir = '/nfs/a65/pmcjs/CMIP6/AerChemMIP/'+model+'/'+pertrun+'/histSST/'
    deltadir = '/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/deltas/'+pertrun+'/piClim-histaer/'
    if has_rfmip:
        deltadir = '/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'_AerChemMIP'+'/deltas/'+pertrun+'/piClim-histaer/'
    mkdir_p(deltadir)
    
    # make deltas and ERF
    rsdt_base = iris.load(basedir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    if make_aprp:
        rlutcs_base = iris.load(basedir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_base)
        rlutcs_base = rlutcs_base.concatenate()[0]
        rsds_base = iris.load(basedir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_base)
        rsds_base = rsds_base.concatenate()[0]
        rsus_base = iris.load(basedir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_base)
        rsus_base = rsus_base.concatenate()[0]
        rsutcs_base = iris.load(basedir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_base)
        rsutcs_base = rsutcs_base.concatenate()[0]
        clt_base = iris.load(basedir+'clt*.nc')
        equalise_attributes(clt_base)
        clt_base = clt_base.concatenate()[0]
    if make_od550aer:
        od550aer_base = iris.load(basedir+'od550aer_*.nc', 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
        equalise_attributes(od550aer_base)
        od550aer_base = od550aer_base.concatenate()[0]

    rsdt_pert = iris.load(pertdir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    if make_aprp:
        rsds_pert = iris.load(pertdir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_pert)
        unify_time_units(rsds_pert)
        rsds_pert = rsds_pert.concatenate()[0]
        rsus_pert = iris.load(pertdir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_pert)
        unify_time_units(rsus_pert)
        rsus_pert = rsus_pert.concatenate()[0]
        rlutcs_pert = iris.load(pertdir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_pert)
        unify_time_units(rlutcs_pert)
        rlutcs_pert = rlutcs_pert.concatenate()[0]
        rsutcs_pert = iris.load(pertdir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_pert)
        unify_time_units(rsutcs_pert)
        rsutcs_pert = rsutcs_pert.concatenate()[0]
        clt_pert = iris.load(pertdir+'clt*.nc')
        equalise_attributes(clt_pert)
        unify_time_units(clt_pert)
        clt_pert = clt_pert.concatenate()[0]
    if make_od550aer:
        od550aer_pert = iris.load(pertdir+'od550aer_*.nc', "atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
        equalise_attributes(od550aer_pert)
        unify_time_units(od550aer_pert)
        od550aer_pert = od550aer_pert.concatenate()[0]

    nmonths = rsdt_pert.data.shape[0]
    nmonthsbase = rsdt_base.data.shape[0]
    nyears = nmonths//12
    nyearsbase = nmonthsbase//12
    print(nmonths)

    erf_sw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_sw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_sw_data = (rsdt_pert.data-rsut_pert.data)-(rsdt_base.data-rsut_base.data)
    erf_lw_data = (-rlut_pert.data)-(-rlut_base.data)

    if make_aprp:
        erf_sw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_lw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_sw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_lw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_sw_cs_data = (rsdt_pert.data-rsutcs_pert.data)-(rsdt_base.data-rsutcs_base.data)
        erf_lw_cs_data = (-rlutcs_pert.data)-(-rlutcs_base.data)

    if make_od550aer:
        od550aer_data = (od550aer_pert.data-od550aer_base.data)

    erf_sw = iris.cube.Cube(
        erf_sw_data,
        var_name = 'erf_sw',
        long_name = "Shortwave effective radiative forcing at TOA",
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_sw.coord('time').attributes = {}
    erf_sw.coord('latitude').attributes = {}
    erf_sw.coord('longitude').attributes = {}
    iris.save(erf_sw, deltadir+'erf_sw.nc')
    
    erf_lw = iris.cube.Cube(
        erf_lw_data,
        var_name = 'erf_lw',
        long_name = 'Longwave effective radiative forcing at TOA',
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_lw.coord('time').attributes = {}
    erf_lw.coord('latitude').attributes = {}
    erf_lw.coord('longitude').attributes = {}
    iris.save(erf_lw, deltadir+'erf_lw.nc')

    erf = erf_lw + erf_sw.data
    erf.var_name = 'erf'
    erf.long_name = 'Effective radiative forcing at TOA'
    iris.save(erf, deltadir+'erf.nc')

    if make_aprp:
        erf_sw_cs = iris.cube.Cube(
            erf_sw_cs_data,
            var_name = 'erf_sw_cs',
            long_name = "Shortwave clear-sky effective radiative forcing at TOA",
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_sw_cs.coord('time').attributes = {}
        erf_sw_cs.coord('latitude').attributes = {}
        erf_sw_cs.coord('longitude').attributes = {}
        iris.save(erf_sw_cs, deltadir+'erf_sw_cs.nc')
        
        erf_lw_cs = iris.cube.Cube(
            erf_lw_cs_data,
            var_name = 'erf_lw_cs',
            long_name = 'Longwave clear-sky effective radiative forcing at TOA',
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_lw_cs.coord('time').attributes = {}
        erf_lw_cs.coord('latitude').attributes = {}
        erf_lw_cs.coord('longitude').attributes = {}
        iris.save(erf_lw_cs, deltadir+'erf_lw_cs.nc')
        
        erf_cs = erf_lw_cs + erf_sw_cs.data
        erf_cs.var_name = 'erf_cs'
        erf_cs.long_name = 'Clear-sky effective radiative forcing at TOA'
        iris.save(erf_cs, deltadir+'erf_cs.nc')
    
    if make_od550aer:
        od550aer = iris.cube.Cube(
            od550aer_data,
            var_name = 'od550aer',
            long_name = 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles',
            units = '1',
            dim_coords_and_dims=[(od550aer_pert.coord('time'), 0), (od550aer_pert.coord('latitude'), 1), (od550aer_pert.coord('longitude'), 2)]
        )
        od550aer.coord('time').attributes = {}
        od550aer.coord('latitude').attributes = {}
        od550aer.coord('longitude').attributes = {}
        iris.save(od550aer, deltadir+'od550aer.nc')

    # easiest to do by working directly with the data
    adjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
    adjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'
    if has_rfmip:
        adjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'_AerChemMIP/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
        adjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'_AerChemMIP/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'

    mkdir_p(adjalldir)
    mkdir_p(adjglobalmeanyearmeandir)

    # global means, statistics, etc
    variables = ['erf_sw', 'erf_lw', 'erf']
    if make_aprp:
        variables.extend(['erf_sw_cs', 'erf_lw_cs', 'erf_cs'])
    if make_od550aer:
        variables.extend(['od550aer'])
    for var in variables:
        cube = iris.load_cube(deltadir+var+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, adjglobalmeanyearmeandir+var+'.nc')
    return


def make_deltas(baserun, pertrun, model, expt, make_aprp, make_od550aer):
    basedir = '/nfs/a65/pmcjs/CMIP6/RFMIP/'+model+'/'+baserun+'/piClim-control/'
    pertdir = '/nfs/a65/pmcjs/CMIP6/RFMIP/'+model+'/'+pertrun+'/piClim-'+expt+'/'
    deltadir = '/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/deltas/'+pertrun+'/piClim-'+expt+'/'
    mkdir_p(deltadir)
    
    # make deltas and ERF
    rsdt_base = iris.load(basedir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    if make_aprp:
        rlutcs_base = iris.load(basedir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_base)
        rlutcs_base = rlutcs_base.concatenate()[0]
        rsds_base = iris.load(basedir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_base)
        rsds_base = rsds_base.concatenate()[0]
        rsus_base = iris.load(basedir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_base)
        rsus_base = rsus_base.concatenate()[0]
        rsutcs_base = iris.load(basedir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_base)
        rsutcs_base = rsutcs_base.concatenate()[0]
        clt_base = iris.load(basedir+'clt*.nc')
        equalise_attributes(clt_base)
        clt_base = clt_base.concatenate()[0]
    if make_od550aer:
        od550aer_base = iris.load(basedir+'od550aer_*.nc', 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
        equalise_attributes(od550aer_base)
        od550aer_base = od550aer_base.concatenate()[0]

    rsdt_pert = iris.load(pertdir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    if make_aprp:
        rsds_pert = iris.load(pertdir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
        equalise_attributes(rsds_pert)
        unify_time_units(rsds_pert)
        rsds_pert = rsds_pert.concatenate()[0]
        rsus_pert = iris.load(pertdir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
        equalise_attributes(rsus_pert)
        unify_time_units(rsus_pert)
        rsus_pert = rsus_pert.concatenate()[0]
        rlutcs_pert = iris.load(pertdir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
        equalise_attributes(rlutcs_pert)
        unify_time_units(rlutcs_pert)
        rlutcs_pert = rlutcs_pert.concatenate()[0]
        rsutcs_pert = iris.load(pertdir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
        equalise_attributes(rsutcs_pert)
        unify_time_units(rsutcs_pert)
        rsutcs_pert = rsutcs_pert.concatenate()[0]
        clt_pert = iris.load(pertdir+'clt*.nc')
        equalise_attributes(clt_pert)
        unify_time_units(clt_pert)
        clt_pert = clt_pert.concatenate()[0]
    if make_od550aer:
        od550aer_pert = iris.load(pertdir+'od550aer_*.nc', "atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
        equalise_attributes(od550aer_pert)
        unify_time_units(od550aer_pert)
        od550aer_pert = od550aer_pert.concatenate()[0]

    nmonths = rsdt_pert.data.shape[0]
    nmonthsbase = rsdt_base.data.shape[0]
    nyears = nmonths//12
    nyearsbase = nmonthsbase//12
    print(nmonths)

    erf_sw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_sw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    erf_lw_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))

    if make_aprp:
        erf_sw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_lw_cs_data = np.zeros((nmonths, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_sw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
        erf_lw_cs_temp = np.zeros((nyearsbase, 12, rsdt_pert.shape[1], rsdt_pert.shape[2]))
    if make_od550aer:
        od550aer_data = np.zeros((nmonths, od550aer_pert.shape[1], od550aer_pert.shape[2]))
        od550aer_temp = np.zeros((nyearsbase, 12, od550aer_pert.shape[1], od550aer_pert.shape[2]))


    for i in range(nyears):
        for j in range(nyearsbase):
            erf_sw_temp[j,:,:,:] = (rsdt_pert.data[i*12:i*12+12,...]-rsut_pert.data[i*12:i*12+12,...]) - (rsdt_base.data[j*12:j*12+12,...]-rsut_base.data[j*12:j*12+12,...])
            erf_lw_temp[j,:,:,:] = (-rlut_pert.data[i*12:i*12+12,...]) - (-rlut_base.data[j*12:j*12+12,...])
            if make_aprp:
                erf_sw_cs_temp[j,:,:,:] = (rsdt_pert.data[i*12:i*12+12,...]-rsutcs_pert.data[i*12:i*12+12,...]) - (rsdt_base.data[j*12:j*12+12,...]-rsutcs_base.data[j*12:j*12+12,...])
                erf_lw_cs_temp[j,:,:,:] = (-rlutcs_pert.data[i*12:i*12+12,...]) - (-rlutcs_base.data[j*12:j*12+12,...])
            if make_od550aer:
                od550aer_temp[j,:,:,:] = od550aer_pert.data[i*12:i*12+12,...] - od550aer_base.data[j*12:j*12+12,...]
        erf_sw_data[i*12:i*12+12,:,:] = np.mean(erf_sw_temp, axis=0)
        erf_lw_data[i*12:i*12+12,:,:] = np.mean(erf_lw_temp, axis=0)
        if make_aprp:
            erf_sw_cs_data[i*12:i*12+12,:,:] = np.mean(erf_sw_cs_temp, axis=0)
            erf_lw_cs_data[i*12:i*12+12,:,:] = np.mean(erf_lw_cs_temp, axis=0)
        if make_od550aer:
            od550aer_data[i*12:i*12+12,:,:] = np.mean(od550aer_temp, axis=0)

    erf_sw = iris.cube.Cube(
        erf_sw_data,
        var_name = 'erf_sw',
        long_name = "Shortwave effective radiative forcing at TOA",
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_sw.coord('time').attributes = {}
    erf_sw.coord('latitude').attributes = {}
    erf_sw.coord('longitude').attributes = {}
    iris.save(erf_sw, deltadir+'erf_sw.nc')
    
    erf_lw = iris.cube.Cube(
        erf_lw_data,
        var_name = 'erf_lw',
        long_name = 'Longwave effective radiative forcing at TOA',
        units = "W m-2",
        dim_coords_and_dims=[(rsdt_pert.coord('time'), 0), (rsdt_pert.coord('latitude'), 1), (rsdt_pert.coord('longitude'), 2)]
    )
    erf_lw.coord('time').attributes = {}
    erf_lw.coord('latitude').attributes = {}
    erf_lw.coord('longitude').attributes = {}
    iris.save(erf_lw, deltadir+'erf_lw.nc')

    erf = erf_lw + erf_sw.data
    erf.var_name = 'erf'
    erf.long_name = 'Effective radiative forcing at TOA'
    iris.save(erf, deltadir+'erf.nc')

    if make_aprp:    
        erf_sw_cs = iris.cube.Cube(
            erf_sw_cs_data,
            var_name = 'erf_sw_cs',
            long_name = "Shortwave clear-sky effective radiative forcing at TOA",
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_sw_cs.coord('time').attributes = {}
        erf_sw_cs.coord('latitude').attributes = {}
        erf_sw_cs.coord('longitude').attributes = {}
        iris.save(erf_sw_cs, deltadir+'erf_sw_cs.nc')
        
        erf_lw_cs = iris.cube.Cube(
            erf_lw_cs_data,
            var_name = 'erf_lw_cs',
            long_name = 'Longwave clear-sky effective radiative forcing at TOA',
            units = "W m-2",
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        erf_lw_cs.coord('time').attributes = {}
        erf_lw_cs.coord('latitude').attributes = {}
        erf_lw_cs.coord('longitude').attributes = {}
        iris.save(erf_lw_cs, deltadir+'erf_lw_cs.nc')
        
        erf_cs = erf_lw_cs + erf_sw_cs.data
        erf_cs.var_name = 'erf_cs'
        erf.long_name = 'Clear-sky effective radiative forcing at TOA'
        iris.save(erf_cs, deltadir+'erf_cs.nc')

    if make_od550aer:
        od550aer = iris.cube.Cube(
            od550aer_data,
            var_name = 'od550aer',
            long_name = 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles',
            units = '1',
            dim_coords_and_dims=[(od550aer_pert.coord('time'), 0), (od550aer_pert.coord('latitude'), 1), (od550aer_pert.coord('longitude'), 2)]
        )
        od550aer.coord('time').attributes = {}
        od550aer.coord('latitude').attributes = {}
        od550aer.coord('longitude').attributes = {}
        iris.save(od550aer, deltadir+'od550aer.nc')

    
    # easiest to do by working directly with the data
    deltadir = '/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/deltas/'+pertrun+'/piClim-'+expt+'/'
    adjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/all/'
    adjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/globalmeanyearmean/'

    mkdir_p(adjalldir)
    mkdir_p(adjglobalmeanyearmeandir)

    # global means, statistics, etc
    variables = ['erf_sw', 'erf_lw', 'erf']
    if make_aprp:
        variables.extend(['erf_sw_cs', 'erf_lw_cs', 'erf_cs'])
    if make_od550aer:
        variables.extend(['od550aer'])
    for var in variables:
        cube = iris.load_cube(deltadir+var+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, adjglobalmeanyearmeandir+var+'.nc')
    return


def make_aprp_rfmip(baserun, pertrun, model, expt):
    basedir = '/nfs/a65/pmcjs/CMIP6/RFMIP/'+model+'/'+baserun+'/piClim-control/'
    pertdir = '/nfs/a65/pmcjs/CMIP6/RFMIP/'+model+'/'+pertrun+'/piClim-'+expt+'/'
    aprpadjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/all/'
    aprpadjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/globalmeanyearmean/'

    mkdir_p(aprpadjalldir)
    mkdir_p(aprpadjglobalmeanyearmeandir)

    clt_base = iris.load(basedir+'clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_base)
    clt_base = clt_base.concatenate()[0]
    rsds_base = iris.load(basedir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_base)
    rsds_base = rsds_base.concatenate()[0]
    rsus_base = iris.load(basedir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_base)
    rsus_base = rsus_base.concatenate()[0]
    rsdscs_base = iris.load(basedir+'rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_base)
    rsdscs_base = rsdscs_base.concatenate()[0]
    rsuscs_base = iris.load(basedir+'rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_base)
    rsuscs_base = rsuscs_base.concatenate()[0]
    rsdt_base = iris.load(basedir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rlutcs_base = iris.load(basedir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_base)
    rlutcs_base = rlutcs_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    rsutcs_base = iris.load(basedir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_base)
    rsutcs_base = rsutcs_base.concatenate()[0]

    clt_pert = iris.load(pertdir+'clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_pert)
    unify_time_units(clt_pert)
    clt_pert = clt_pert.concatenate()[0]
    rsds_pert = iris.load(pertdir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_pert)
    unify_time_units(rsds_pert)
    rsds_pert = rsds_pert.concatenate()[0]
    rsus_pert = iris.load(pertdir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_pert)
    unify_time_units(rsus_pert)
    rsus_pert = rsus_pert.concatenate()[0]
    rsdscs_pert = iris.load(pertdir+'rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_pert)
    unify_time_units(rsdscs_pert)
    rsdscs_pert = rsdscs_pert.concatenate()[0]
    rsuscs_pert = iris.load(pertdir+'rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_pert)
    unify_time_units(rsuscs_pert)
    rsuscs_pert = rsuscs_pert.concatenate()[0]
    rsdt_pert = iris.load(pertdir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rlutcs_pert = iris.load(pertdir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_pert)
    unify_time_units(rlutcs_pert)
    rlutcs_pert = rlutcs_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    rsutcs_pert = iris.load(pertdir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_pert)
    unify_time_units(rsutcs_pert)
    rsutcs_pert = rsutcs_pert.concatenate()[0]

    nrpt = int(np.ceil(rsds_pert.data.shape[0] / rsds_base.data.shape[0]))
    nmonths = rsds_pert.data.shape[0]
    nmonthsbase = rsds_base.data.shape[0]
    nyears = nmonths//12
    nyearsbase = nmonthsbase//12
    print(nrpt, nmonths)
    shortlist = ['ERFariSW', 'ERFaciSW', 'albedo', 'ERFariLW', 'ERFaciLW', 't9']
    aprp = {}
    aprp_temp = {}

    for item in shortlist:
        aprp[item] = np.zeros((nmonths, rsds_pert.data.shape[1], rsds_pert.data.shape[2]))

    for i in range(nyears):
        for item in shortlist:
            aprp_temp[item] = np.zeros((nyearsbase, 12, rsds_base.data.shape[1], rsds_base.data.shape[2]))
        for j in range(nyearsbase):
            base = {
                'rsdt'   : rsdt_base.data[j*12:j*12+12,...],
                'rsus'   : rsus_base.data[j*12:j*12+12,...],
                'rsds'   : rsds_base.data[j*12:j*12+12,...],
                'clt'    : clt_base.data[j*12:j*12+12,...],
                'rsdscs' : rsdscs_base.data[j*12:j*12+12,...],
                'rsuscs' : rsuscs_base.data[j*12:j*12+12,...],
                'rsut'   : rsut_base.data[j*12:j*12+12,...],
                'rsutcs' : rsutcs_base.data[j*12:j*12+12,...],
                'rlut'   : rlut_base.data[j*12:j*12+12,...],
                'rlutcs' : rlutcs_base.data[j*12:j*12+12,...],
            }
            pert = {
                'rsdt'   : rsdt_pert.data[i*12:i*12+12,...],
                'rsus'   : rsus_pert.data[i*12:i*12+12,...],
                'rsds'   : rsds_pert.data[i*12:i*12+12,...],
                'clt'    : clt_pert.data[i*12:i*12+12,...],
                'rsdscs' : rsdscs_pert.data[i*12:i*12+12,...],
                'rsuscs' : rsuscs_pert.data[i*12:i*12+12,...],
                'rsut'   : rsut_pert.data[i*12:i*12+12,...],
                'rsutcs' : rsutcs_pert.data[i*12:i*12+12,...],
                'rlut'   : rlut_pert.data[i*12:i*12+12,...],
                'rlutcs' : rlutcs_pert.data[i*12:i*12+12,...]
            }
        
            #for key in base.keys():
            #    print (key, base[key].shape, pert[key].shape)
        
            aprp_output = calc_aprp(base, pert, lw=True)
            component_longnames={
                't1'      : 'Change in shortwave clear-sky radiation with respect to surface albedo',
                't2'      : 'Change in shortwave clear-sky radiation with respect to aerosol scattering',
                't3'      : 'Change in shortwave clear-sky radiation with respect to aerosol absorption',
                't4'      : 'Change in shortwave overcast radiation with respect to surface albedo',
                't5'      : 'Change in shortwave overcast radiation with respect to aerosol scattering',
                't6'      : 'Change in shortwave overcast radiation with respect to aerosol absorption',
                't7'      : 'Change in shortwave overcast radiation with respect to cloud scattering',
                't8'      : 'Change in shortwave overcast radiation with respect to cloud absorption',
                't9'      : 'Change in shortwave radiation with respect to cloud fraction',
                'ERFariSW': 'Shortwave effective radiative forcing due to aerosol-radiation interactions',
                'ERFaciSW': 'Shortwave effective radiative forcing due to aerosol-cloud interactions',
                'albedo'  : 'Shortwave effective radiative forcing due to surface albedo',
                'ERFariLW': 'Longwave effective radiative forcing due to aerosol-radiation interactions',
                'ERFaciLW': 'Longwave effective radiative forcing due to aerosol-cloud interactions',
                't2_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol scattering assuming clear sky',
                't3_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol absorption assuming clear sky',
                'ERFariSWclr': 'Shortwave effective radiative forcing due to aerosol-radiation interactions assuming clear sky'
            }

            for item in shortlist:
                aprp_temp[item][j,:,:,:] = aprp_output[item]
        for item in shortlist:
            aprp[item][12*i:12*i+12,:,:] = np.mean(aprp_temp[item],axis=0)

    for component in aprp.keys():
        cube = iris.cube.Cube(
            aprp[component],
            var_name = component,
            long_name = component_longnames[component],
            units = 'W m-2',
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        iris.save(cube, aprpadjalldir+'aprp_'+component+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, aprpadjglobalmeanyearmeandir+'aprp_'+component+'.nc')
    return

def make_aprp_e3sm(baserun, pertrun, model, expt):
    basedir = '/nfs/a65/pmcjs/E3SM/%s/Amon/' % baserun
    pertdir = '/nfs/a65/pmcjs/E3SM/%s/Amon/' % pertrun
    aprpadjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
    aprpadjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'

    mkdir_p(aprpadjalldir)
    mkdir_p(aprpadjglobalmeanyearmeandir)

    clt_base = iris.load(basedir+'clt/gr/v20200528/clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_base)
    clt_base = clt_base.concatenate()[0]
    rsds_base = iris.load(basedir+'rsds/gr/v20200528/rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_base)
    rsds_base = rsds_base.concatenate()[0]
    rsus_base = iris.load(basedir+'rsus/gr/v20200528/rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_base)
    rsus_base = rsus_base.concatenate()[0]
    rsdscs_base = iris.load(basedir+'rsdscs/gr/v20200528/rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_base)
    rsdscs_base = rsdscs_base.concatenate()[0]
    rsuscs_base = iris.load(basedir+'rsuscs/gr/v20200528/rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_base)
    rsuscs_base = rsuscs_base.concatenate()[0]
    rsdt_base = iris.load(basedir+'rsdt/gr/v20200528/rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut/gr/v20200528/rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rlutcs_base = iris.load(basedir+'rlutcs/gr/v20200528/rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_base)
    rlutcs_base = rlutcs_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut/gr/v20200528/rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    rsutcs_base = iris.load(basedir+'rsutcs/gr/v20200528/rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_base)
    rsutcs_base = rsutcs_base.concatenate()[0]

    clt_pert = iris.load(pertdir+'clt/gr/v20200528/clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_pert)
    unify_time_units(clt_pert)
    clt_pert = clt_pert.concatenate()[0]
    rsds_pert = iris.load(pertdir+'rsds/gr/v20200528/rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_pert)
    unify_time_units(rsds_pert)
    rsds_pert = rsds_pert.concatenate()[0]
    rsus_pert = iris.load(pertdir+'rsus/gr/v20200528/rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_pert)
    unify_time_units(rsus_pert)
    rsus_pert = rsus_pert.concatenate()[0]
    rsdscs_pert = iris.load(pertdir+'rsdscs/gr/v20200528/rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_pert)
    unify_time_units(rsdscs_pert)
    rsdscs_pert = rsdscs_pert.concatenate()[0]
    rsuscs_pert = iris.load(pertdir+'rsuscs/gr/v20200528/rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_pert)
    unify_time_units(rsuscs_pert)
    rsuscs_pert = rsuscs_pert.concatenate()[0]
    rsdt_pert = iris.load(pertdir+'rsdt/gr/v20200528/rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut/gr/v20200528/rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rlutcs_pert = iris.load(pertdir+'rlutcs/gr/v20200528/rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_pert)
    unify_time_units(rlutcs_pert)
    rlutcs_pert = rlutcs_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut/gr/v20200528/rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    rsutcs_pert = iris.load(pertdir+'rsutcs/gr/v20200528/rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_pert)
    unify_time_units(rsutcs_pert)
    rsutcs_pert = rsutcs_pert.concatenate()[0]

    nmonths = rsds_pert.data.shape[0]
    nyears = nmonths//12
    print(nmonths)
#    shortlist = ['ERFariSW', 'ERFaciSW', 'albedo', 'ERFariLW', 'ERFaciLW', 't9']
    aprp = {}

#    for item in shortlist:
#        aprp[item] = np.zeros((nmonths, rsds_pert.data.shape[1], rsds_pert.data.shape[2]))

    base = {
        'rsdt'   : rsdt_base.data,
        'rsus'   : rsus_base.data,
        'rsds'   : rsds_base.data,
        'clt'    : clt_base.data,
        'rsdscs' : rsdscs_base.data,
        'rsuscs' : rsuscs_base.data,
        'rsut'   : rsut_base.data,
        'rsutcs' : rsutcs_base.data,
        'rlut'   : rlut_base.data,
        'rlutcs' : rlutcs_base.data,
    }
    pert = {
        'rsdt'   : rsdt_pert.data,
        'rsus'   : rsus_pert.data,
        'rsds'   : rsds_pert.data,
        'clt'    : clt_pert.data,
        'rsdscs' : rsdscs_pert.data,
        'rsuscs' : rsuscs_pert.data,
        'rsut'   : rsut_pert.data,
        'rsutcs' : rsutcs_pert.data,
        'rlut'   : rlut_pert.data,
        'rlutcs' : rlutcs_pert.data
    }
    
    aprp = calc_aprp(base, pert, lw=True)
    component_longnames={
        't1'      : 'Change in shortwave clear-sky radiation with respect to surface albedo',
        't2'      : 'Change in shortwave clear-sky radiation with respect to aerosol scattering',
        't3'      : 'Change in shortwave clear-sky radiation with respect to aerosol absorption',
        't4'      : 'Change in shortwave overcast radiation with respect to surface albedo',
        't5'      : 'Change in shortwave overcast radiation with respect to aerosol scattering',
        't6'      : 'Change in shortwave overcast radiation with respect to aerosol absorption',
        't7'      : 'Change in shortwave overcast radiation with respect to cloud scattering',
        't8'      : 'Change in shortwave overcast radiation with respect to cloud absorption',
        't9'      : 'Change in shortwave radiation with respect to cloud fraction',
        'ERFariSW': 'Shortwave effective radiative forcing due to aerosol-radiation interactions',
        'ERFaciSW': 'Shortwave effective radiative forcing due to aerosol-cloud interactions',
        'albedo'  : 'Shortwave effective radiative forcing due to surface albedo',
        'ERFariLW': 'Longwave effective radiative forcing due to aerosol-radiation interactions',
        'ERFaciLW': 'Longwave effective radiative forcing due to aerosol-cloud interactions',
        't2_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol scattering assuming clear sky',
        't3_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol absorption assuming clear sky',
        'ERFariSWclr': 'Shortwave effective radiative forcing due to aerosol-radiation interactions assuming clear sky'
    }

    for component in aprp.keys():
        cube = iris.cube.Cube(
            aprp[component],
            var_name = component,
            long_name = component_longnames[component],
            units = 'W m-2',
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        iris.save(cube, aprpadjalldir+'aprp_'+component+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, aprpadjglobalmeanyearmeandir+'aprp_'+component+'.nc')
    return

def make_aprp_aerchemmip(baserun, pertrun, model, expt, has_rfmip):
    basedir = '/nfs/a65/pmcjs/CMIP6/AerChemMIP/'+model+'/'+baserun+'/histSST-piAer/'
    pertdir = '/nfs/a65/pmcjs/CMIP6/AerChemMIP/'+model+'/'+pertrun+'/histSST/'
    aprpadjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
    aprpadjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'
    if has_rfmip:
        aprpadjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'_AerChemMIP/adjust/aprp/'+pertrun+'/piClim-histaer/all/'
        aprpadjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'_AerChemMIP/adjust/aprp/'+pertrun+'/piClim-histaer/globalmeanyearmean/'
    mkdir_p(aprpadjalldir)
    mkdir_p(aprpadjglobalmeanyearmeandir)

    clt_base = iris.load(basedir+'clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_base)
    clt_base = clt_base.concatenate()[0]
    rsds_base = iris.load(basedir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_base)
    rsds_base = rsds_base.concatenate()[0]
    rsus_base = iris.load(basedir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_base)
    rsus_base = rsus_base.concatenate()[0]
    rsdscs_base = iris.load(basedir+'rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_base)
    rsdscs_base = rsdscs_base.concatenate()[0]
    rsuscs_base = iris.load(basedir+'rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_base)
    rsuscs_base = rsuscs_base.concatenate()[0]
    rsdt_base = iris.load(basedir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rlutcs_base = iris.load(basedir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_base)
    rlutcs_base = rlutcs_base.concatenate()[0]
    rsut_base = iris.load(basedir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    rsutcs_base = iris.load(basedir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_base)
    rsutcs_base = rsutcs_base.concatenate()[0]

    clt_pert = iris.load(pertdir+'clt_*.nc', 'cloud_area_fraction')
    equalise_attributes(clt_pert)
    unify_time_units(clt_pert)
    clt_pert = clt_pert.concatenate()[0]
    rsds_pert = iris.load(pertdir+'rsds_*.nc', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_pert)
    unify_time_units(rsds_pert)
    rsds_pert = rsds_pert.concatenate()[0]
    rsus_pert = iris.load(pertdir+'rsus_*.nc', 'surface_upwelling_shortwave_flux_in_air')
    equalise_attributes(rsus_pert)
    unify_time_units(rsus_pert)
    rsus_pert = rsus_pert.concatenate()[0]
    rsdscs_pert = iris.load(pertdir+'rsdscs_*.nc', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_pert)
    unify_time_units(rsdscs_pert)
    rsdscs_pert = rsdscs_pert.concatenate()[0]
    rsuscs_pert = iris.load(pertdir+'rsuscs_*.nc', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_pert)
    unify_time_units(rsuscs_pert)
    rsuscs_pert = rsuscs_pert.concatenate()[0]
    rsdt_pert = iris.load(pertdir+'rsdt_*.nc', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    unify_time_units(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'rlut_*.nc', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    unify_time_units(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rlutcs_pert = iris.load(pertdir+'rlutcs_*.nc', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_pert)
    unify_time_units(rlutcs_pert)
    rlutcs_pert = rlutcs_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'rsut_*.nc', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    unify_time_units(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    rsutcs_pert = iris.load(pertdir+'rsutcs_*.nc', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_pert)
    unify_time_units(rsutcs_pert)
    rsutcs_pert = rsutcs_pert.concatenate()[0]

    nmonths = rsds_pert.data.shape[0]
    nyears = nmonths//12
    print(nmonths)
    aprp = {}

    base = {
        'rsdt'   : rsdt_base.data,
        'rsus'   : rsus_base.data,
        'rsds'   : rsds_base.data,
        'clt'    : clt_base.data,
        'rsdscs' : rsdscs_base.data,
        'rsuscs' : rsuscs_base.data,
        'rsut'   : rsut_base.data,
        'rsutcs' : rsutcs_base.data,
        'rlut'   : rlut_base.data,
        'rlutcs' : rlutcs_base.data,
    }
    pert = {
        'rsdt'   : rsdt_pert.data,
        'rsus'   : rsus_pert.data,
        'rsds'   : rsds_pert.data,
        'clt'    : clt_pert.data,
        'rsdscs' : rsdscs_pert.data,
        'rsuscs' : rsuscs_pert.data,
        'rsut'   : rsut_pert.data,
        'rsutcs' : rsutcs_pert.data,
        'rlut'   : rlut_pert.data,
        'rlutcs' : rlutcs_pert.data
    }
    
    aprp = calc_aprp(base, pert, lw=True)
    component_longnames={
        't1'      : 'Change in shortwave clear-sky radiation with respect to surface albedo',
        't2'      : 'Change in shortwave clear-sky radiation with respect to aerosol scattering',
        't3'      : 'Change in shortwave clear-sky radiation with respect to aerosol absorption',
        't4'      : 'Change in shortwave overcast radiation with respect to surface albedo',
        't5'      : 'Change in shortwave overcast radiation with respect to aerosol scattering',
        't6'      : 'Change in shortwave overcast radiation with respect to aerosol absorption',
        't7'      : 'Change in shortwave overcast radiation with respect to cloud scattering',
        't8'      : 'Change in shortwave overcast radiation with respect to cloud absorption',
        't9'      : 'Change in shortwave radiation with respect to cloud fraction',
        'ERFariSW': 'Shortwave effective radiative forcing due to aerosol-radiation interactions',
        'ERFaciSW': 'Shortwave effective radiative forcing due to aerosol-cloud interactions',
        'albedo'  : 'Shortwave effective radiative forcing due to surface albedo',
        'ERFariLW': 'Longwave effective radiative forcing due to aerosol-radiation interactions',
        'ERFaciLW': 'Longwave effective radiative forcing due to aerosol-cloud interactions',
        't2_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol scattering assuming clear sky',
        't3_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol absorption assuming clear sky',
        'ERFariSWclr': 'Shortwave effective radiative forcing due to aerosol-radiation interactions assuming clear sky'
    }

    for component in aprp.keys():
        cube = iris.cube.Cube(
            aprp[component],
            var_name = component,
            long_name = component_longnames[component],
            units = 'W m-2',
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        iris.save(cube, aprpadjalldir+'aprp_'+component+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, aprpadjglobalmeanyearmeandir+'aprp_'+component+'.nc')
    return

def make_aprp_from_pp(baserun, pertrun, model, expt):
    basedir = '/nfs/a65/pmcjs/RFMIP/tim/piClim-control-notESGF/'
    pertdir = '/nfs/a65/pmcjs/CMIP6/RFMIP/'+model+'/'+pertrun+'/piClim-'+expt+'/'
#    pertdir = '/nfs/a65/pmcjs/RFMIP/tim/piClim-histaer-notESGF/'
    aprpadjalldir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/all/'
    aprpadjglobalmeanyearmeandir='/nfs/a65/pmcjs/RFMIP/RFMIP-ERF/'+model+'/adjust/aprp/'+pertrun+'/piClim-'+expt+'/globalmeanyearmean/'

    mkdir_p(aprpadjalldir)
    mkdir_p(aprpadjglobalmeanyearmeandir)

    clt_base = iris.load(basedir+'*.pp', 'cloud_area_fraction')
    equalise_attributes(clt_base)
    clt_base = clt_base.concatenate()[0]
    rsds_base = iris.load(basedir+'*.pp', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_base)
    rsds_base = rsds_base.concatenate()[0]
    rsns_base = iris.load(basedir+'*.pp', 'surface_net_downward_shortwave_flux')
    equalise_attributes(rsns_base)
    rsns_base = rsns_base.concatenate()[0]
    rsdscs_base = iris.load(basedir+'*.pp', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_base)
    rsdscs_base = rsdscs_base.concatenate()[0]
    rsuscs_base = iris.load(basedir+'*.pp', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_base)
    rsuscs_base = rsuscs_base.concatenate()[0]
    rsdt_base = iris.load(basedir+'*.pp', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_base)
    rsdt_base = rsdt_base.concatenate()[0]
    rlut_base = iris.load(basedir+'*.pp', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_base)
    rlut_base = rlut_base.concatenate()[0]
    rlutcs_base = iris.load(basedir+'*.pp', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_base)
    rlutcs_base = rlutcs_base.concatenate()[0]
    rsut_base = iris.load(basedir+'*.pp', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_base)
    rsut_base = rsut_base.concatenate()[0]
    rsutcs_base = iris.load(basedir+'*.pp', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_base)
    rsutcs_base = rsutcs_base.concatenate()[0]

    rsus_base = rsds_base - rsns_base

    clt_pert = iris.load(pertdir+'*.pp', 'cloud_area_fraction')
    equalise_attributes(clt_pert)
    clt_pert = clt_pert.concatenate()[0]
    rsds_pert = iris.load(pertdir+'*.pp', 'surface_downwelling_shortwave_flux_in_air')
    equalise_attributes(rsds_pert)
    rsds_pert = rsds_pert.concatenate()[0]
    rsns_pert = iris.load(pertdir+'*.pp', 'surface_net_downward_shortwave_flux')
    equalise_attributes(rsns_pert)
    rsns_pert = rsns_pert.concatenate()[0]
    rsdscs_pert = iris.load(pertdir+'*.pp', 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsdscs_pert)
    rsdscs_pert = rsdscs_pert.concatenate()[0]
    rsuscs_pert = iris.load(pertdir+'*.pp', 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky')
    equalise_attributes(rsuscs_pert)
    rsuscs_pert = rsuscs_pert.concatenate()[0]
    rsdt_pert = iris.load(pertdir+'*.pp', 'toa_incoming_shortwave_flux')
    equalise_attributes(rsdt_pert)
    rsdt_pert = rsdt_pert.concatenate()[0]
    rlut_pert = iris.load(pertdir+'*.pp', 'toa_outgoing_longwave_flux')
    equalise_attributes(rlut_pert)
    rlut_pert = rlut_pert.concatenate()[0]
    rlutcs_pert = iris.load(pertdir+'*.pp', 'toa_outgoing_longwave_flux_assuming_clear_sky')
    equalise_attributes(rlutcs_pert)
    rlutcs_pert = rlutcs_pert.concatenate()[0]
    rsut_pert = iris.load(pertdir+'*.pp', 'toa_outgoing_shortwave_flux')
    equalise_attributes(rsut_pert)
    rsut_pert = rsut_pert.concatenate()[0]
    rsutcs_pert = iris.load(pertdir+'*.pp', 'toa_outgoing_shortwave_flux_assuming_clear_sky')
    equalise_attributes(rsutcs_pert)
    rsutcs_pert = rsutcs_pert.concatenate()[0]

    rsus_pert = rsds_pert - rsns_pert

    nrpt = int(np.ceil(rsds_pert.data.shape[0] / rsds_base.data.shape[0]))
    nmonths = rsds_pert.data.shape[0]
    nmonthsbase = rsds_base.data.shape[0]
    nyears = nmonths//12
    nyearsbase = nmonthsbase//12
    print(nrpt, nmonths)
    shortlist = ['ERFariSW', 'ERFaciSW', 'albedo', 'ERFariLW', 'ERFaciLW', 't9']
    aprp = {}
    aprp_temp = {}

    for item in shortlist:
        aprp[item] = np.zeros((nmonths, rsds_pert.data.shape[1], rsds_pert.data.shape[2]))

    for i in range(nyears):
        for item in shortlist:
            aprp_temp[item] = np.zeros((nyearsbase, 12, rsds_base.data.shape[1], rsds_base.data.shape[2]))
        for j in range(nyearsbase):
            base = {
                'rsdt'   : rsdt_base.data[j*12:j*12+12,...],
                'rsus'   : rsus_base.data[j*12:j*12+12,...],
                'rsds'   : rsds_base.data[j*12:j*12+12,...],
                'clt'    : clt_base.data[j*12:j*12+12,...]*100,
                'rsdscs' : rsdscs_base.data[j*12:j*12+12,...],
                'rsuscs' : rsuscs_base.data[j*12:j*12+12,...],
                'rsut'   : rsut_base.data[j*12:j*12+12,...],
                'rsutcs' : rsutcs_base.data[j*12:j*12+12,...],
                'rlut'   : rlut_base.data[j*12:j*12+12,...],
                'rlutcs' : rlutcs_base.data[j*12:j*12+12,...],
            }
            pert = {
                'rsdt'   : rsdt_pert.data[i*12:i*12+12,...],
                'rsus'   : rsus_pert.data[i*12:i*12+12,...],
                'rsds'   : rsds_pert.data[i*12:i*12+12,...],
                'clt'    : clt_pert.data[i*12:i*12+12,...]*100,
                'rsdscs' : rsdscs_pert.data[i*12:i*12+12,...],
                'rsuscs' : rsuscs_pert.data[i*12:i*12+12,...],
                'rsut'   : rsut_pert.data[i*12:i*12+12,...],
                'rsutcs' : rsutcs_pert.data[i*12:i*12+12,...],
                'rlut'   : rlut_pert.data[i*12:i*12+12,...],
                'rlutcs' : rlutcs_pert.data[i*12:i*12+12,...]
            }
        
            #for key in base.keys():
            #    print (key, base[key].shape, pert[key].shape)
        
            aprp_output = calc_aprp(base, pert, lw=True)
            component_longnames={
                't1'      : 'Change in shortwave clear-sky radiation with respect to surface albedo',
                't2'      : 'Change in shortwave clear-sky radiation with respect to aerosol scattering',
                't3'      : 'Change in shortwave clear-sky radiation with respect to aerosol absorption',
                't4'      : 'Change in shortwave overcast radiation with respect to surface albedo',
                't5'      : 'Change in shortwave overcast radiation with respect to aerosol scattering',
                't6'      : 'Change in shortwave overcast radiation with respect to aerosol absorption',
                't7'      : 'Change in shortwave overcast radiation with respect to cloud scattering',
                't8'      : 'Change in shortwave overcast radiation with respect to cloud absorption',
                't9'      : 'Change in shortwave radiation with respect to cloud fraction',
                'ERFariSW': 'Shortwave effective radiative forcing due to aerosol-radiation interactions',
                'ERFaciSW': 'Shortwave effective radiative forcing due to aerosol-cloud interactions',
                'albedo'  : 'Shortwave effective radiative forcing due to surface albedo',
                'ERFariLW': 'Longwave effective radiative forcing due to aerosol-radiation interactions',
                'ERFaciLW': 'Longwave effective radiative forcing due to aerosol-cloud interactions',
                't2_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol scattering assuming clear sky',
                't3_clr'  : 'Change in shortwave clear-sky radiation with respect to aerosol absorption assuming clear sky',
                'ERFariSWclr': 'Shortwave effective radiative forcing due to aerosol-radiation interactions assuming clear sky'
            }

            for item in shortlist:
                aprp_temp[item][j,:,:,:] = aprp_output[item]
        for item in shortlist:
            aprp[item][12*i:12*i+12,:,:] = np.mean(aprp_temp[item],axis=0)

    for component in aprp.keys():
        cube = iris.cube.Cube(
            aprp[component],
            var_name = component,
            long_name = component_longnames[component],
            units = 'W m-2',
            dim_coords_and_dims=[(rsds_pert.coord('time'), 0), (rsds_pert.coord('latitude'), 1), (rsds_pert.coord('longitude'), 2)]
        )
        iris.save(cube, aprpadjalldir+'aprp_'+component+'.nc')
        iris.coord_categorisation.add_year(cube, 'time')
        cube_year = cube.aggregated_by('year', iris.analysis.MEAN)
        if not cube_year.coord('latitude').has_bounds():
            cube_year.coord('latitude').guess_bounds()
        if not cube_year.coord('longitude').has_bounds():
            cube_year.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube_year)
        cube_gmym = cube_year.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
        iris.save(cube_gmym, aprpadjglobalmeanyearmeandir+'aprp_'+component+'.nc')
    return


def driver(baserun, pertrun, model, expt, from_pp=False, make_aprp=True, make_od550aer=False, e3sm=False, aerchemmip=False, has_rfmip=False):
    if from_pp:
        make_deltas_from_pp(baserun, pertrun, model, expt)
    elif e3sm:
        make_deltas_e3sm(baserun, pertrun, model, expt, make_aprp)
    elif aerchemmip:
        make_deltas_aerchemmip(baserun, pertrun, model, expt, make_aprp, make_od550aer, has_rfmip)
    else:
        make_deltas(baserun, pertrun, model, expt, make_aprp, make_od550aer)
    if make_aprp:
        if from_pp:
            make_aprp_from_pp(baserun, pertrun, model, expt)
        elif e3sm:
            make_aprp_e3sm(baserun, pertrun, model, expt)
        elif aerchemmip:
            make_aprp_aerchemmip(baserun, pertrun, model, expt, has_rfmip)
        else:
            make_aprp_rfmip(baserun, pertrun, model, expt)
    return
 
