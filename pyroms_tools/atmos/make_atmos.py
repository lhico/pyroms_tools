import xarray as xr
from pyroms_tools.utils  import atmos_forcing as af 
from pyroms_tools import utils as ut
import numpy as np
import pandas as pd
from netCDF4 import date2num
from scipy import ndimage
import os
import glob
import yaml
import argparse
from datetime import datetime as dtt


"""
This script preprocesses the era5 winds for roms

required variables for the script (single levels product of ERA5):

ATTENTION TO THE UNITS! You can switch fluxes for mean fluxes provided you adjust the scale factors
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
varinfo.yml in ROMS

"10m_u_component_of_wind",     # used in dQdSST
"10m_v_component_of_wind",     # used in dQdSST
"2m_temperature",              # used in dQdSST/air_density
"sea_level_pressure",          # used in dQdSST/air_density
"sea_surface_temperature",     # used in masks
"surface_latent_heat_flux",    # used in heat_flux
"surface_sensible_heat_flux",  # used in heat_flux
"surface_net_solar_radiation",   # used in heat_flux
"surface_net_thermal_radiation", # used in heat_flux
"total_precipitation",         # used in net_freshwater
"evaporation",                 # used in net_freshwater

required variables for the script (pressure levels product of ERA5)?
"specific_humidity"            # used in air_density

"""

def mixing_ratio(specific_humidity):
    q = specific_humidity
    return q / (1 - q)

def virtual_temperature(temperature, mixing_ratio):
    tv = temperature * (1 + 0.608 * mixing_ratio)
    return tv

def air_density(pressure, virtual_temperature):
    rho = pressure / (287 * virtual_temperature)
    return rho

def mask_sst(ds):
    """
    Create a mask for sea surface temperature.

    Parameters:
    ds (xarray.Dataset): The dataset containing the sea surface temperature.
    dicts (dict): Dictionary containing the keys and variable ids
                    of the sea surface temperature.
                    The keys are:
                        - 'sea_surface_temperature'
    """
    mask = ds['sea_surface_temperature']
    condition = ~np.isnan(mask.isel(time=0))
    return xr.where(condition, 0, 1)

def net_freshwater(ds, scale=1):
    """
    Calculate the net freshwater flux from various components in the dataset.

    Parameters:
    ds (xarray.Dataset): The dataset containing freshwater flux components.
    dicts (dict): Dictionary containing the keys and variable ids
                    of the freshwater flux components.  
                    The keys are:
                        - 'evaporation'
                        - 'total_precipitation' 
    """
    evaporation = ds['evaporation']
    precipitation = ds['total_precipitation']
    return (evaporation + precipitation) * scale

def heat_flux(ds, scale=1):
    """
    Calculate the total heat flux from various components in the dataset.

    Parameters:
    ds (xarray.Dataset): The dataset containing heat flux components.
    dicts (dict): Dictionary containing the keys and variable ids 
                  of the heat flux components.
                  The keys are:
                    - 'surface_sensible_heat_flux'
                    - 'surface_latent_heat_flux'
                    - 'surface_net_solar_radiation'
                    - 'surface_net_thermal_radiation'
    scale (float): Scaling factor for the total heat flux (default 1).

    Returns:
    xarray.DataArray: The total heat flux.
    """
    qs = ds['surface_sensible_heat_flux']
    ql = ds['surface_latent_heat_flux']
    snsr = ds['surface_net_solar_radiation']
    sntr = ds['surface_net_thermal_radiation']

    return (qs + ql + snsr + sntr) * scale   


def select_time_range(nc, timei, timef):
    tstart_num = date2num(dtt.strptime(timei, '%Y-%m-%dT%H:%M:%S'), nc.time.units)
    tfinal_num = date2num(dtt.strptime(timef, '%Y-%m-%dT%H:%M:%S'), nc.time.units)
    nc = nc.sel(time=slice(tstart_num, tfinal_num))
    return nc



def compute_dQSST(nc):
    u10 = nc['10m_u_component_of_wind']
    v10 = nc['10m_v_component_of_wind']
    d2m = nc['2m_dewpoint_temperature']
    sp  = nc['sea_level_pressure']
    t2m = nc['2m_temperature']

    wspd = (u10 ** 2 + v10 ** 2) ** 0.5
    w = mixing_ratio(d2m)
    tv = virtual_temperature(t2m, w)
    rhoa = air_density(sp, tv)
    dqdsst = af.dQdT(wspd, sp, t2m, rhoair=rhoa)
    nc['dQdSST'] = dqdsst.dQdT()#.isel(level=0)
    return nc

def process_variable(varname, ncout, tref1, metadata,mask):
    rename = metadata[varname]['outputName']
    var = ncout[[varname]]
    var = af.extrapolating_era5(ncout, varname, None, extrapolate_method='laplace', dst=ncout, mask=mask)
    
    if varname == 'sp':
        aux = af.extrapolating_era5(ncout, varname, None, extrapolate_method='xesmf', dst=ncout, mask=abs(ncout['mask'] - 2 + 1))
        var.sp.values = aux.sp.values
        var['sp'] = var['sp'].fillna(var.sp.mean().values)
        aux = ndimage.gaussian_filter(var.sp.values, 1)
        var.sp.values = aux

    var = var.reindex(lat=list(reversed(var.lat.values)))
    var = var.rename_vars({varname: rename,})
    # var = var.rename_dims({'latitude': 'lat', 'longitude': 'lon'})
    var[rename].attrs.update({
        'coordinates': 'lon lat',
        'units': metadata[varname]['units'],
        'scale_factor': 1,
        'add_offset': 0
    })

    var = var.assign_coords({'time': tref1})
    var['time'].attrs['units'] = 'days since 1990-01-01 00:00:00'
    var.attrs['time'] = metadata[varname]['time']
    var.attrs['coordinates'] = 'lon lat'

    if rename in ['SST', 'temp']:
        var[rename].values -= 273.15

    print(rename)
    var.attrs = {}
    return var

def get_standard_dict():
    """
    Convert a dictionary to a list of standard names. These are variables
    from the source dataset that will be used as forcing in ROMS.

    Parameters:
    dicts (dict): Dictionary containing the keys and variable ids.

    Returns:
    list: List of standard names.
    """
    
    s_dicts = {
      '10m_u_component_of_wind': 'u10',       # used in dQdSST
      '10m_v_component_of_wind': 'v10',       # used in dQdSST
      '2m_temperature': 't2m',                # used in dQdSST/air_density
      'sea_level_pressure': 'sp',            # used in dQdSST/air_density
      'turbulent_east_surface_stress': 'metss',
    #   'turbulent_north_surface_stress': 'mntss',
    }
    
    return s_dicts


def switched_dict_items_keys(dicts):
    return {value: key for key, value in dicts.items()}

def main():
	#if __name__ == '__main__':
    """
    Bulk fluxes used in ROMS are momentum, heat and freshwater fluxes.
    
    Momentum fluxes are either calculated from the 10m wind components or passed
    directly from the source dataset. Heat fluxes are calculated from the surface
    sensible and latent heat fluxes, and the surface net solar and thermal radiation. 
    Freshwater fluxes are calculated from the evaporation and precipitation.
    
    The heat flux sensitivity to sea surface temperature is also calculated.
    
    Since the source dataset may be relatively low resolution and can bias the ocean
    simulation, the data is extrapolated to the coast using the Laplace method our
    a simple nearest neighbor method.

    The original source dataset is the ERA5 reanalysis from ECMWF.
    """
    parser = argparse.ArgumentParser(description="Process ERA5 data for ROMS")
    parser.add_argument('--config', type=str, help='Path to the YAML configuration file')
    args = parser.parse_args()
    config_path = args.config



    # Load YAML configuration
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Access the parameters from the YAML
    era_config = config['default']['era']
    fpath_pressure = era_config['press']
    fpath_single = era_config['single']
    odir = era_config['outputdir']
    timei, timef, dt = era_config['date']

    # load the dictionary with the renaming of the variables
    # we opted to have the full names of each field for clarity
    era_rename_vars = config['default']['era']['fields']
    era_rename_coords = config['default']['era']['coordinates']

    # Load datasets
    fpath = glob.glob(fpath_single)
    print(fpath)
    nc = [xr.open_dataset(f, decode_times=False) for f in fpath]
    # nc2 = xr.open_mfdataset(fpath[5:9],combine='by_coords', compat='override', decode_times=False)
    # nc3 = xr.open_mfdataset(fpath[10:],combine='by_coords', compat='override', decode_times=False)
    nc = xr.merge(nc)
    print(nc)
    # standardizing coordinate names to time, lat, lon
    aux = switched_dict_items_keys(era_rename_coords)
    nc = nc.rename(aux)

    # standardizing variable names
    aux = switched_dict_items_keys(era_rename_vars) # convert to long name
#     standard_dicts = standard_dict()
#     aux1 = [s_dicts]
    nc = nc.rename(aux)

#     # Select time range
    nc= select_time_range(nc, timei, timef)

#     # Create time reference

    tref = pd.date_range(start=timei, periods=nc.time.size, freq='1H')
    tref1 = date2num(tref.to_pydatetime(), 'days since 1990-01-01 00:00:00')



    # Compute fluxes
    nc['shflux'] = heat_flux(nc, scale=1 / 3600) # heat flux
    nc['swflux'] = -net_freshwater(nc, scale=1e-3 / 3600)

    nc['mask'] = mask_sst(nc) #create a land/water mask
    # Compute sensitivity
    nc = compute_dQSST(nc)

    # Get the standard dictionary (convert nc to have short variable names)
    standard_dict = get_standard_dict()
    ncout = nc.rename(standard_dict)
    standard_dict_switched = switched_dict_items_keys(standard_dict)
    #creating a list that will be used in the data processing
    standard_list = [f for f in standard_dict_switched]
    # standard_list.append('mask')
    standard_list.append('dQdSST')
    standard_list.append('shflux')
    standard_list.append('swflux')



    # Process each variable and save to netCDF
    # WARNING: notice that af.variables_list contain the actual
    # naming conventions that will be used in the final output
    for varname in standard_list:
        if glob.glob(f'test_{varname}_{tref.year[0]}-{tref.month[0]:02d}.nc'):
            print(f'{varname} already saved')
            continue
        var = process_variable(varname, ncout, tref1, af.variables_list, nc['mask'])
        var.to_netcdf(f'test_{varname}_{tref.year[0]}-{tref.month[0]:02d}.nc', format='NETCDF4')

    # Save the final output
    dsout = xr.open_mfdataset(f'test_*_{tref.year[0]}-{tref.month[0]:02d}.nc', compat='override')
    dsout.drop_vars('step').to_netcdf(odir)
    # dsout
    # os.system(f'rm test_*_{tref.year[0]}-{tref.month[0]:02d}.nc')

if __name__ == "__main__":


     main()
