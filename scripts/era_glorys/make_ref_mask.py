import xarray as xr
from utils import configs
import os.path as osp
import matplotlib.pyplot as plt

datapath = '../../data/roms_files'
fpath = '../../data/aux_data/GLO-MFC_001_030_mask_bathy.nc'

nclocal = xr.open_dataset(osp.join(datapath, 'sbb_grid_roms.nc'))

extent = configs.area_params['pcse1']
# extent = configs.area_params['wSAtl0']
# fpath = '/home/otel/Desktop/GLO-MFC_001_030_mask_bathy.nc'
nc = xr.open_dataset(fpath)
ncut = nc.sel(longitude=slice(extent['lonW'], extent['lonE']),
              latitude=slice(extent['latS'], extent['latN']))
ncut.to_netcdf(osp.join(datapath, 'glo_mask_bathy2.nc'))

# use the following figures to write xrange and yrange in the xrange and yrange
# fields in areaparams in config.py (xrange and yrange are the index that
# limits the grid area)
plt.figure(1)
plt.contourf(ncut.longitude, ncut.latitude, ncut.deptho)
plt.contourf(nclocal.lon_rho-360, nclocal.lat_rho, nclocal.h)
