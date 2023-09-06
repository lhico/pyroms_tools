import xarray as xr
from utils.atmos_forcing import dQdT, extrapolating_era5, variables_list
from utils import utils as ut
import numpy as np
import pandas as pd
from matplotlib import rcParams
from netCDF4 import date2num, num2date
from scipy import ndimage
import sys
import os

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
    return q/(1-q)

def virtual_temperature(temperature, mixing_ratio):
    tv = temperature*(1 + 0.608*mixing_ratio)
    return tv

def air_density(pressure, virtual_temperature):
    rho = pressure/ (287*virtual_temperature)
    return rho

def mask_sst(ds, varb='sst'):
    condition = ~np.isnan(ds[varb].isel(time=0))
    mask = xr.where(condition, 0, 1)
    return mask

def net_freshwater(ds,
                   evaporation='e',
                   precipitation='tp',
                   scale=1,
                   ):

    return (ds[evaporation] + ds[precipitation])*scale

def heat_flux(ds,
              qs='sshf', # surface sensible heat flux
              ql='slhf', # surface latent heat flux
              ssr='ssr', # surface short radiation 
              str='str', # surface thermal radiation
              scale=1):
    return (ds[qs] + ds[ql] + ds[ssr] + ds[str])*scale   



# -- gets  the information from the config file -- #
# getting the referemce domain from shell 
if len(sys.argv) > 1:
    reference = sys.argv[1]
else:
    reference = 'ceresIV_2.012'

dicts = ut._get_dict_paths(f'{os.path.dirname(__file__)}/../configs/grid_config_esmf.txt')[reference]


# dicts resolver
fpath_single   = dicts['frc.era.single']
fpath_pressure = dicts['frc.era.press']
time0          = dicts['frc.date'][0]
timef          = dicts['frc.date'][1]
dt             = dicts['frc.date'][2]
odir           = dicts['frc.outputdir']


# month = 5
nc = xr.open_dataset(fpath_single, decode_times=False)
nc1 = xr.open_dataset(fpath_pressure, decode_times=False)
nc['q'] = nc1['q']
sstK = True


time0 = num2date(nc.time[0].values, nc.time.attrs['units'])
tref = pd.date_range(start=str(time0), periods=nc.time.size, freq='1H')
tref1 = date2num(tref.to_pydatetime(), 'days since 1990-01-01 00:00:00')



metadata = variables_list


nc['mask'] = mask_sst(nc, varb='sst')

# -- heat flux -- #
nc['shflux'] = heat_flux(nc,
                        qs='sshf', # surface sensible heat flux
                        ql='slhf', # surface latent heat flux
                        ssr='ssr', # surface short radiation 
                        str='str', # surface thermal radiation)
                        scale=1/3600)

# -- fresh water fluxes -- #
nc['swflux'] = net_freshwater(nc,
                            evaporation='e',
                            precipitation='tp',
                            scale=1e-3)

# -- surface net heat flux sensitivity to SST -- #
wspd = (nc['u10']**2 + nc['u10']**2)**0.5
w    = mixing_ratio(nc['q'])
tv   = virtual_temperature(nc['t2m'], w)
rhoa = air_density(nc['sp'], tv)

dqdsst = dQdT(wspd,nc['sp'], nc['t2m'], rhoair=rhoa)
nc['dQdSST']= dqdsst.dQdT().isel(level=0)

ncout = nc[['mask', 'shflux', 'swflux', 'metss', 'mntss', 'dQdSST', 'sst','sp']]

for varname in ['sp','sst', 'metss', 'mntss', 'shflux', 'swflux', 'dQdSST']:
    rename = metadata[varname]['outputName']
    var = extrapolating_era5(ncout, varname, None, 
                                extrapolate_method='laplace', dst=ncout, mask=nc['mask'])

    if varname=='sp':
        aux = extrapolating_era5(ncout, varname, None, 
                                extrapolate_method='xesmf', dst=ncout, mask=abs(nc['mask']-2+1))
        var.sp.values = aux.sp.values
        var['sp'] = var['sp'].fillna(var.sp.mean().values)
        aux = ndimage.gaussian_filter(var.sp.values,1)
        var.sp.values = aux
    var = var.reindex(latitude=list(reversed(var.latitude.values)))
    var = var.rename({varname: rename,
                    'latitude': 'lat',
                    'longitude': 'lon'})
    var[rename].attrs['coordinates'] = 'lon lat' 
    var[rename].attrs['units'] = metadata[varname]['units']
    var[rename].attrs['scale_factor'] = 1
    var[rename].attrs['add_offset'] = 0



    # var = var.assign_coords({'time': tref})
    # # var['time'].attrs['units'] = 'days since 1970-01-01'

    var = var.assign_coords({'time': tref1})
    var['time'].attrs['units'] = 'days since 1990-01-01 00:00:00'

    var.attrs['time'] = metadata[varname]['time']
    var.attrs['coordinates'] = 'lon lat'

    if rename=='SST':
        var[rename].values -= 273.15
    # var = var.assign_coords(time=pd.date_range(start=ncout.time.values[0], freq='1H', periods=ncout.time.size))
    print(rename)
    var.attrs = {}
    var.to_netcdf(f'test_{rename}_{tref.year[0]}-{tref.month[0]:02d}.nc', format='NETCDF4')

dsout = xr.open_mfdataset(f'test_*_{tref.year[0]}-{tref.month[0]:02d}.nc')
dsout.to_netcdf(f'{odir}')
os.system(f'rm test_*_{tref.year[0]}-{tref.month[0]:02d}.nc')
