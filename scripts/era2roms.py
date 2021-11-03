import xarray as xr
import os.path as osp
import glob
import datetime as dtt
from utils.barnier95 import dQdT
import utils as ut
import os
os.environ["PYROMS_GRIDID_FILE"] = "/home/lhico/pyroms_tools/configs/gridid.txt"

def scaling(x):
    sf = x.attrs['scale_factor']
    ao = x.attrs['add_offset']
    return sf*x + ao



# -- gets  the information from the config file -- #
reference = 'pbs_202109_glorys'
dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

rename       = dicts['atmos.rename']
atmos_single = dicts['atmos.single']
atmos_pressr = dicts['atmos.pressure']
atmostr     = dicts['atmos.str']
outdir       = dicts['output_dir']


for ia, ip in zip(atmos_single, atmos_pressr):
    tstr = osp.basename(ia)
    tstr = str(dtt.datetime.strptime(tstr, atmostr))
    tstr = tstr.replace(' ', 'T')
    # -- read netcdfs -- #
    nc = xr.open_dataset(osp.join(ia), decode_cf=False)
    nc1 = xr.open_dataset(osp.join(ip), decode_cf=False)

    # -- append humidity to the main netcdf file -- #
    nc['q'] = nc1['q']

    # -- renaming netcdf variables to roms requirments -- #
    keylist = list(rename.values())[:-2]
    nc = nc.assign(longitude=lambda x: x.longitude + 360)
    ncrename = nc.rename(rename)

    # -- adjusting latitude to be monotonic increasing (roms doesn't like)
    # decreasing values --#
    ncrename['lat'] = ncrename['lat'][::-1]
    ncrename = ncrename.assign_coords({'time': ncrename['time'].values/24})
    ncrename.time.attrs['units'] = 'days since 1900-01-01 00:00:00.0'

    ncrename['Pair'].attrs['scale_factor'] = 0.01 * ncrename['Pair'].attrs['scale_factor']
    ncrename['Pair'].attrs['add_offset'] = 0.01 * ncrename['Pair'].attrs['add_offset']
    ncrename['Pair'].attrs['units'] = 'mbar'

    ncrename['Tair'].attrs['units'] = 'degC'
    ncrename['Tair'].attrs['add_offset'] = ncrename['Tair'].attrs['add_offset']-273.15

    ncrename['SST'].attrs['units'] = 'degC'
    ncrename['SST'].attrs['add_offset'] = ncrename['SST'].attrs['add_offset']-273.15

    ncrename['Qair'].attrs['scale_factor'] = 1e3 * ncrename['Qair'].attrs['scale_factor']
    ncrename['Qair'].attrs['units'] = 'g kg**-1'




    for strvarb in keylist:
        ncrename[strvarb].attrs['coordinates'] = 'lon lat'
        ncrename[strvarb].attrs['time'] = 'time'
        # adjusting fields to change of latitude coordinates #
        ncrename[strvarb].values = ncrename[strvarb][:, ::-1, :].values

        wind = (scaling(ncrename['Uwind'])**2 + scaling(ncrename['Vwind'])**2)**0.5
        p = scaling(ncrename['Pair'])
        T = scaling(ncrename['Tair']) + 273.15
        rhoair = 1.16

        q = dQdT(wind, p, T)
        dQdSST = q.dQdT()
    ncrename['dQdSST'] = dQdSST
    ncrename['dQdSST'].attrs['coordinates'] = 'lon lat'
    ncrename['dQdSST'].attrs['time'] = 'time'
    ncrename['dQdSST'].attrs['scale_factor'] = 1
    ncrename['dQdSST'].attrs['add_offset'] = 0

    ncrename.to_netcdf(osp.join(outdir, f'era5_{tstr}.nc'))
