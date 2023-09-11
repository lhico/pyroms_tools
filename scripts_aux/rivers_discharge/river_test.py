import pandas as pd
import xarray as xr
from dry_package.plot_schemes import maps
import cartopy.crs as ccrs
import numpy as np
"""
script (test) to produce river output in roms for cananeia estuary
"""


def grdc_laplata():
    # read GRDC ata
    ds1 = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/data/rivers/GRDC/laplata_RS_SC_RJ/2023-05-02_19-37/GRDC-Daily.nc')

    # read birocchi data
    ds2 = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/data/rivers/birocchi_etal_2023/river_valogrande_measured_model.nc')
    ds2 = ds2.assign_coords(id=[1])

    # Adjusting birocchi data to have the same length as GRDC
    # Create a new DataArray with NaN values for the extended period
    new_time = pd.date_range(ds1.time[0].values, ds1.time[-1].values, freq='D')
    new_data = np.empty(len(new_time))
    new_data[:] = np.nan
    new_da = xr.DataArray(new_data.reshape(new_data.size,1), coords=[new_time, ds2.coords['id']], dims=['time', 'id'])
    extended_da = xr.concat([ds2['model'], new_da], dim='time')
    extended_da = extended_da.sortby('time')
    extended_da = extended_da.drop_duplicates(dim='time')
    extended_da = extended_da.to_dataset()
    extended_da = extended_da.rename({'model':'runoff_mean'})

    # we will output an xarray dataset and we are preparing the dictionaries here
    v = 'runoff_mean'
    rivers = {
        'laplata'    : (['time'], ds1[v].sel(id=3265601).values +  # parana
                                  ds1[v].sel(id=3469050).values,  # uruguai
                       # ds1[v].sel(id=3469910),# negro
                       ),
        'lagoapatos' : (['time'], ds1[v].sel(id=3653400).values +  # jacui
                                  ds1[v].sel(id=3653770).values),  # gavatai
        'tubarao'    : (['time'], ds1[v].sel(id=3653450).values,),  
        'itajai'     : (['time'], ds1[v].sel(id=3653352).values,),  
        'itapocu'    : (['time'], ds1[v].sel(id=3653250).values,),
        'mambucaba'  : (['time'], ds1[v].sel(id=3653985).values),
        'iguape'     : (['time'], extended_da[v].values.squeeze(),)
    }

    dims = {'time': ds1.time.values}

    ds = xr.Dataset(data_vars=rivers, coords=dims)
    return ds


def fill_with_climatology(dataarray):
    q = dataarray.copy(deep=True)  # avoid pointers issue
    qclim = q.groupby('time.dayofyear').mean() # calculate daily average

    # identify nans that will be replaced with avg values
    inan  = np.where(np.isnan(q))[0]
    
    # replace nans at given daysofyear with daily avg values
    for i in range(366):
        idayofyear = q['time.dayofyear'] == i+1
        isnanq     = np.isnan(q)
        q.values[(idayofyear) & (isnanq)] = qclim.values[i]
    return q


def roms_river(roms, dslist, templist, saltlist, tstart, tend, dirlist=None, signlist=None):
    roms = roms.copy()
    roms = roms.sel(river_time=slice(tstart,tend))

    
    for i, ds1 in enumerate(dslist):
        dsaux = ds1.sel(time=slice(tstart,tend))
        dsaux = dsaux.resample(time=timefreq).mean()
        roms.river_transport.values[:, i]= dsaux.values

    for i,t in enumerate(templist):
        roms.river_temp.values[:,:, i] = t

    for i,s in enumerate(saltlist):
        roms.river_salt.values[:,:, i] = s

    # roms = roms.drop(['river_direction'])
    if dirlist is not None:
        roms['river_direction'] = ('river', dirlist)
    else:
        roms['river_direction'] = roms.river_direction.astype(int)
   
    if signlist is not None:
        roms.river_sign.values[:] = signlist

    roms.river_direction.values = roms.river_direction.astype(float)

    return roms



# read data
roms_river_template  = '/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools_new/data/treated_deep4/river_template.nc'
# roms_river_template2  = '/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools_new/data/treated_deep4/nested/river_template.nc'

grid = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools_new/data/treated/bacia_santos_final.nc')

roms  = xr.open_dataset(roms_river_template)
# roms2 = xr.open_dataset(roms_river_template2)


# # if needed, plot to check points
ax = maps.make_overall(extent=[-50,-40,-30,-20])
ax.coastlines()
# x = grid.lon_rho[roms2.river_Eposition[0].astype(int), roms2.river_Xposition[0].astype(int)]
# y = grid.lat_rho[roms2.river_Eposition[0].astype(int), roms2.river_Xposition[0].astype(int)]
# ax.scatter(roms2.lon, roms2.lat)
# ax.scatter(x,y)



ds1 = grdc_laplata()
ds2 = ds1.copy()
for varname in ds1.data_vars:
    ds2[varname] = fill_with_climatology(ds1[varname])



tstart = roms.river_time.min()
tend   = roms.river_time.max()
timefreq = '1D'

dslist = [ds2['laplata']/5, ds2['laplata']/5,ds2['laplata']/5, ds2['laplata']/5,ds2['laplata']/5, ds2['lagoapatos']/3, ds2['lagoapatos']/3, ds2['lagoapatos']/3]
roms = roms_river(roms,
                  dslist,
                  [17,17,17,17,17,20,20],
                  [0,0,0],
                  tstart,tend)

# dslist = [ds2['tubarao'], ds2['itajai'], ds2['itapocu'], ds2['iguape'], ds2['mambucaba']]
# roms2 = roms_river(roms2,
#                    dslist,
#                    np.ones(5)*20,
#                    np.zeros(5)*0,
#                    tstart,tend,
#                    dirlist=np.zeros(5),
#                    signlist=np.zeros(5))

roms.to_netcdf('rivers.nc')
# roms2.to_netcdf('rivers_nested.nc')