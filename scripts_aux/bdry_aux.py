import xarray as xr
import numpy as np
import datetime as dtt
import pandas as pd

t = pd.date_range(start='2019-08-01T12:00:00', end='2019-09-01T12:00:00', freq='1D')
nc = xr.open_dataset('pbs_202109_smooth_bdry_2019-08-01T12:00:00.nc')
for i in t[1:]:
     a = dtt.datetime.strftime(i, 'pbs_202109_smooth_bdry_%Y-%m-%dT12:00:00.nc')
     nc1 = nc.copy()
     nc1 = nc1.assign(ocean_time=[i])
     nc1.to_netcdf(a)