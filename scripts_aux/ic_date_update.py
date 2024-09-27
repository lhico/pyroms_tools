import xarray as xr
import pandas as pd
import sys

nc = xr.open_dataset('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_ic_2.012.nc')
nc = nc.assign_coords(ocean_time=pd.date_range(start='2010-05-08T12', periods=1, freq='1D'))
nc.to_netcdf('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_ic_2.012_update.nc')