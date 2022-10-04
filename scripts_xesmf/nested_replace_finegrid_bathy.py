import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.backend_bases import MouseButton
import xesmf as xe
import os.path as osp

# replaces the bathymetry from the ncchange with a high resolution bathymetry
masklist = ['h', 'mask_u', 'mask_v']
lonlist  = ['lon_rho', 'lon_u', 'lon_v']
latlist  = ['lat_rho', 'lat_u', 'lat_v']

# open datasets
ncref    = xr.open_dataset('/home/otel/Desktop/teste/pyroms_tools/data/raw/gebco.nc')
ncchange = xr.open_dataset('/home/otel/Desktop/teste/pyroms_tools/data/treated/roms_grid_smooth_swatl_2022_nested.nc')
outname = osp.basename('/home/otel/Desktop/teste/pyroms_tools/data/treated/roms_grid_smooth_swatl_2022_nested.nc')

# renaming variables to xesmf standards
ncchange = ncchange.rename(lon_rho='lon', lat_rho='lat')
ncchange = ncchange.set_coords(['lon', 'lat'])


# subsampling gebco - full resolution may be to heavy for memory
# TODO: this could be problematic, since the subsampling at certain 
# points may not be represent a given location properly
ncref = ncref.sel(lon=slice(ncchange.lon.min(),ncchange.lon.max(),4),
                  lat=slice(ncchange.lat.min(),ncchange.lat.max(),4),)


# interpolated coarse grid onto fine
regridder = xe.Regridder(ncref.elevation, ncchange.h, 'nearest_s2d', extrap_method='nearest_s2d')
interpolated = regridder(ncref)  # interpolating

# replacing values
interpolated.elevation.values[interpolated.elevation.values>=0] = 5
interpolated.elevation.values[interpolated.elevation.values<=-2500] = -2500
interpolated.elevation.values[interpolated.elevation.values<0] *= -1  # roms uses positive bathymetry

ncchange['h'].values = interpolated.elevation
ncchange['hraw'].values = [interpolated.elevation]

# renaming coordinate names back to roms conventions
ncchange = ncchange.rename(lon='lon_rho', lat='lat_rho')

# saving the netcdf with a higher resolution bathymetry
ncchange.to_netcdf(f'{outname[:-3]}_corrected.nc')