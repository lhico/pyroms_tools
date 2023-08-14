import xarray as xr
import pandas as pd

ds = xr.open_mfdataset('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_bdry_2.012*00*nc')
ds = ds.assign_coords(ocean_time=pd.date_range(start='2010-05-01T12', periods=ds.ocean_time.size, freq='1D'))
ds['salt_east'] = ds.salt_east.fillna(34.66913988)
ds['temp_east']=ds.temp_east.fillna(-0.16397166)
ds['u_east'] = ds['u_east'].fillna(0)
ds['v_east'] = ds['v_east'].fillna(0)
ds.to_netcdf('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_bdrys_2.012_2010-05.nc')