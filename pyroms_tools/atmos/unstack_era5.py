import xarray as xr
import glob
import datetime as dtt
fpath = glob.glob('*.grib')


fpath = [f for f in fpath if 'era5' in f]

for f in fpath: 

    with xr.open_dataset(f, engine='cfgrib') as ds:

        for vars in ds.data_vars.keys():
            if 'step' in ds[vars].dims:
                ds1 = ds.rename(time='time0')
                ds2 = ds1[vars].stack(time=('time0', 'step')).reset_index('time')
                ds2 = ds2.assign_coords(time=ds2.valid_time)
                ds3 = ds2.drop_vars(['time0','step','valid_time'])
            else:
                ds2 = ds[vars]
                ds3 = ds2.drop_vars(['valid_time'])
            ds4 = ds3.transpose('time', 'latitude', 'longitude')
            ds5 = ds4.sel(time=slice(dtt.datetime(1993,1,1), dtt.datetime(2020,1,1)))
            ds5.to_netcdf(f.split('.grb')[0]+'.nc')
            ds5.to_dataset().close()
