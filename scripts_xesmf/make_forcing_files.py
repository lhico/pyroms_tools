# code based on d_ecmwf2roms.m from Hernan Arango and John Wilkin

# making sure we are in the right folder to load local scripts
import os
os.chdir("/home/danilo/phd/thesis/chapter3_driving/scripts/pyroms_tools/scripts")

import sys
import xarray as xr
from glob import glob

from utils.atmos_forcing import relative_humidity, dQdT, scaling, extrapolating_era5
from utils.variables_dict_ecmwf import variables_list
from utils import utils as ut

# -- gets  the information from the config file -- #
# getting the referemce domain from shell 
if len(sys.argv) > 1:
    reference = sys.argv[1]
else:
    reference = 'bulkfluxes'

dicts = ut._get_dict_paths('../../../config/grid_config_pyroms.txt')[reference]

preffix = reference
mean_rate = dicts['atmos.meanrate']

# if BULK_FLUXES is activated in your experiment. This will control which variables
# will be processed
bulkfluxes = dicts['atmos.bulkfluxes']

# activate this key if you desire to extrapolate oceanic data over nearshore
# gaps that has land data
extrap = dicts['atmos.extrap']

# -- source folder where ECMWF files are -- #
# get the folder with forcings files downloaded by download_ecmwf.py script
# note that we complement the path with the suffix used to save the ECMWF files in the line #158
experiment_suffix = "" #reference+
datadir = dicts['atmos.listfiles']

# where final forcing files must be saved
savedir = dicts['atmos.savedir']

## -- don't change anything from here, unless you know what you doing -- ##

# variables of interest, based on whether or not BULK_FLUXES is activated
print(f'[SETTINGS] Settings variables list to create forcing files:')
print(f"[SETTINGS] Extrapolation set as {extrap}\n[SETTINGS] BULK FLUXES set as {bulkfluxes}\n[SETTINGS] Mean rate set as {mean_rate}")

if bulkfluxes == True:
    # vars = ['msl', 'u10', 'v10', 't2m', 'tcc', 'sshf', 
    #         'slhf', 'tp', 'str', 'strd', 'ssr', 'q', 
    #         'sst', 'dQdSST']
    vars = ['msl', 'u10', 'v10', 't2m', 'tcc', 'tp', 'ssr', 'q', 'str']
elif bulkfluxes == False:
    if mean_rate:
        vars = ['shflux', 'swflux', 'metss', 'mntss', 'msnswrf', 'dQdSST', 'sst']
    else:
        vars = ['shflux', 'swflux', 'ewss', 'nsss', 'ssr', 'dQdSST', 'sst']
else:
    vars = ['msl', 'u10', 'v10', 't2m', 'tcc', 'sshf', 
            'slhf', 'tp', 'str', 'strd', 'ssr', 'q', 
            'sst', 'dQdSST',
            'shflux', 'swflux', 'ewss', 'nsss', 'ssr']    

print("[LOADING]  source files from ERA5")
# listing all variables' files into a single xarray.Dataset
nfiles = glob(f"{datadir}/*{experiment_suffix}.nc") #{experiment_suffix}/
# read all files into one Dataset 
ds = xr.open_mfdataset(nfiles, decode_cf=False)
sst = xr.open_dataset([f for f in nfiles if 'sea_surface_temperature' in f][0], decode_cf=False)

if extrap == 'xesmf':
    print("[LOADING] destiny grid for xESMF interpolation and extrapolation")
    # destiny grid
    # TODO: adaptar o caminho para o caminho da grade do projeto, segundo config.yml
    dst = xr.open_dataset("/home/danilo/phd/thesis/chapter3_driving/projects/swa/input/grid_SWA5km_1200m.nc")

    # masking properly to xESMF
    dst['mask'] = dst['mask_rho'].copy()
    mask = xr.where(ds['sst'].isel(time=0) == ds['sst'].attrs['missing_value'], 0, 1)
else:
    print("[SETTINGS] mask and dst grid as None")
    mask = None
    dst = None

# start to process
for varname in vars:

    # now begin to process each variable at once, saving into multiple files for 
    # each variable (e.g., sbb_Pair_AUGSEP2019.nc)
    if varname in list(ds.data_vars):
        print("\n --------------------------------")
        print(f"[Processing] converting {varname}")

        # getting metadata
        metadata = variables_list[varname]

        if varname == 't2m':
            print("Converting 2-m air temperature from K to Celsius...")
            ds[varname] = scaling(ds[varname]) - 273.15
            ds[varname].attrs['units'] = 'degC'
            ds[varname].attrs['scale_factor'] = 1
            ds[varname].attrs['add_offset'] = 0
        elif varname == 'sst':
            print("Converting sea surface temperature from K to Celsius...")
            ds[varname] = scaling(ds[varname]) - 273.15
            ds[varname].attrs['units'] = 'degC'
            ds[varname].attrs['scale_factor'] = 1
            ds[varname].attrs['add_offset'] = 0
        else:
            print(f"Scaling {varname}...")
            # just apply the scale factor (conversion of units)
            ds[varname] = scaling(ds[varname]) * metadata['scale']
    else:
        print("\n --------------------------------")
        print(f"[Processing] computing {varname}")
        # computing dQdSST
        if varname == 'dQdSST':
            metadata = variables_list[varname]

            wind = (scaling(ds['u10'])**2 + scaling(ds['v10'])**2) ** 0.5   # scalar wind speed at anemometer level [m s-1]
            p = scaling(ds['msl']) * variables_list['msl']['scale']         # mean sea level pressure [mb]
            #T = scaling(ds['t2m']) + 273.15
            T = scaling(ds['sst'])                                         # sea surface temperature [K]


            dQdSST = dQdT(wind, p, T, rhoair=1.16)
            ds['dQdSST'] = dQdSST.dQdT()
            ds['dQdSST'].attrs['scale_factor'] = 1
            ds['dQdSST'].attrs['add_offset'] = 0

        # computing the net freshwater
        elif varname == 'swflux':
            metadata = variables_list[varname]
            
            if mean_rate:
                # mean rate given in kg m^2 s^-1, so no conversion is needed
                e = -scaling(ds['mer']) #* variables_list['mer']['scale']
                tp= scaling(ds['mtpr']) #* variables_list['mtpr']['scale']
                # (E - TP)/rho_water will give us a swflux in m/s
                ds['swflux'] = (e - tp) * metadata['scale'] # in m s-1
            else:
                # e and tp given as m, so we must convert into kg m^2 s^-1
                e = -scaling(ds['e']) * variables_list['e']['scale']
                tp= scaling(ds['tp']) * variables_list['tp']['scale']
                ds['swflux'] = (e - tp) * metadata['scale'] # in m s-1

           
            ds['swflux'].attrs['scale_factor'] = 1
            ds['swflux'].attrs['add_offset'] = 0
            ds['swflux'].attrs['units'] = metadata['units']
        # computing the net heat flux
        elif varname == 'shflux':
            metadata = variables_list[varname]

            if mean_rate:
                ds['shflux'] = (scaling(ds['msshf']) + 
                                scaling(ds['mslhf']) + 
                                scaling(ds['msnswrf']) + 
                                scaling(ds['msnlwrf'])
                            )
            else:
                ds['shflux'] = (scaling(ds['sshf']) + 
                                scaling(ds['slhf']) + 
                                scaling(ds['str']) + 
                                scaling(ds['ssr'])
                            ) * metadata['scale']

            ds['shflux'].attrs['scale_factor'] = 1
            ds['shflux'].attrs['add_offset'] = 0
            ds['shflux'].attrs['units'] = metadata['units']
            
        # computing the relative humidity
        elif varname == 'q':
            metadata = variables_list[varname]

            # relative humidity requires temperature in Celsius degree
            if ds['t2m'].attrs['units'] != 'degC':
                t2m_C = scaling(ds['t2m']) - 273.15
            else:
                t2m_C = ds['t2m']
            
            d2m_C = scaling(ds['d2m']) - 273.15

            # replacing dewpoint by relative humidity and renaming to
            # the proper ROMS name
            ds['q'] = relative_humidity(t2m_C, d2m_C)
            ds['q'].attrs['scale_factor'] = metadata['scale']

        else:
            print(f'Variable not found in the source file: {varname}')

    # printing maximum and minimum so we can track if something is wrong with the file
    print(f"max = {ds[varname].max().values.item()}\nmin = {ds[varname].min().values.item()}")

    # select variable to save as a separate file
    dsout = ds[varname].to_dataset()

    if extrap:
        print(f"Extrapolating missing values near the coast [method: {extrap}]...")
        field = extrapolating_era5(ds, varname, sst['sst'], 
                                    extrapolate_method=extrap, dst=dst, mask=mask)
        
        # reinserting missing attributes
        field[varname].attrs = ds[varname].attrs
        
        # replacing dsout with the extrapolated field
        dsout[varname] = field[varname]

    # the next commented lines are usefull when using xESMF to interpolate from ERA5 domain to ROMS domains. However, we are
    # not doing this, so I left this piece of code commented just to remember how to implement this functionality.
    # if extrap == 'xesmf':
    #     # rename to ROMS name
    #     dsout = dsout.rename({varname: metadata['outputName'],
    #                             'latitude': 'lat',
    #                             'longitude': 'lon'}
    #                             )

    #     # fix time, latitude, attributes both local and global, and save. Only applicable if we are 
    #     # using the original ECMWF coordinates. If extrapolation was used, then the coordinates are
    #     # the ROMS grid and none inversion is needed.
    #     dsout = dsout.reindex(lat=list(reversed(dsout.lat.values)))
    #     dsout[metadata['outputName']].attrs['coordinates'] = 'lon lat' 

    # else:
    #     dsout = dsout.rename({varname: metadata['outputName']})
    #     dsout[metadata['outputName']].attrs['coordinates'] = 'xi_rho eta_rho'
    
    dsout = dsout.rename({varname: metadata['outputName'],
                            'latitude': 'lat',
                            'longitude': 'lon'}
                            )
    dsout = dsout.reindex(lat=list(reversed(dsout.lat.values)))
    dsout[metadata['outputName']].attrs['coordinates'] = 'lon lat' 
    dsout[metadata['outputName']].attrs['units'] = metadata['units']

    dsout = dsout.assign_coords({'time': dsout['time'].values/24})
    dsout['time'].attrs['units'] = 'days since 1900-01-01 00:00:00.0'
    # replacing temporal axis name
    dsout[metadata['outputName']].attrs['time'] = metadata['time']

    fout = f"{savedir}/{preffix}_{metadata['outputName']}{experiment_suffix}.nc"
    dsout.to_netcdf(fout, format='NETCDF4')
