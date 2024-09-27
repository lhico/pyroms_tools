import pandas as pd
import xarray as xr
import xarray_decorators as xrd
from dry_package.plot_schemes import maps
import cartopy.crs as ccrs
import numpy as np

def grdc_prepare():
    # organize a dataset with the time series of rivers
    ds1 = xr.open_dataset('/home/informatica/Documents/data/rivers/GRDC/laplata_RS_SC_RJ/2023-05-02_19-37/GRDC-Daily.nc')
    ds2 = xr.open_dataset('/home/informatica/Documents/data/rivers/GRDC/paraiba_do_sul/GRDC-Daily.nc')
    ds3 = xr.open_dataset('/home/informatica/Documents/data/rivers/birocchi_etal_2023/river_valogrande_measured_model.nc')
    ds3 = ds3.assign_coords(id=[1])

    ds = xr.merge([ds1,ds2])



    # Adjusting birocchi data to have the same length as GRDC
    # Create a new DataArray with NaN values for the extended period
    new_time = pd.date_range(ds.time[0].values, ds.time[-1].values, freq='D')
    new_data = np.empty(len(new_time))
    new_data[:] = np.nan
    new_da = xr.DataArray(new_data.reshape(new_data.size,1), coords=[new_time, ds3.coords['id']], dims=['time', 'id'])
    extended_da = xr.concat([ds3['model'], new_da], dim='time')
    extended_da = extended_da.sortby('time')
    extended_da = extended_da.drop_duplicates(dim='time')
    extended_da = extended_da.to_dataset()
    extended_da = extended_da.rename({'model':'runoff_mean'})



    # we will output an xarray dataset and we are preparing the dictionaries here
    v = 'runoff_mean'
    rivers = {
        'laplata'    : (['time'], ds[v].sel(id=3265601).values +  # parana
                                  ds[v].sel(id=3469050).values,  # uruguai
                       # ds[v].sel(id=3469910),# negro
                       ),
        'lagoapatos' : (['time'], ds[v].sel(id=3653400).values +  # jacui
                                  ds[v].sel(id=3653770).values),  # gavatai
        'tubarao'    : (['time'], ds[v].sel(id=3653450).values,),  
        'itajai'     : (['time'], ds[v].sel(id=3653352).values,),  
        'itapocu'    : (['time'], ds[v].sel(id=3653250).values,),
        'mambucaba'  : (['time'], ds[v].sel(id=3653985).values),
        'iguape'     : (['time'], extended_da[v].values.squeeze()),
        'macacu'     : (['time'], ds[v].sel(id=3653955).values),
        'macabu'     : (['time'], ds[v].sel(id=3653910).values),
        'paraibasul' : (['time'], ds[v].sel(id=3652890).values), 
    }


    dims = {'time': ds1.time.values}

    ds = xr.Dataset(data_vars=rivers, coords=dims)
    return ds


def roms_river(roms, dslist, templist, saltlist,
               tstart, tend,timefreq, dirlist=None, signlist=None):
    """
    Configure the river file for a ROMS (Regional Ocean Modeling System) simulation.

    Parameters:
    -----------
    roms : xarray.Dataset
        Dataset of a roms river file

    dslist : list of xarray.DataArray
        List containing time series of river discharge for each river.

    templist : list of float
        List of temperatures (constant over time) for each river.

    saltlist : list of float
        List of salinities (constant over time) for each river.

    tstart : numpy.datetime64
        Initial time used to slice the river file.

    tend : numpy.datetime64
        Final time used to slice the river file.

    timefreq : str
        Time frequency for the river data (e.g., '1H' for hourly, '1D' for daily).

    dirlist : list of float, optional
        List of river directions for each river. Default is None.

    signlist : list of float, optional
        List of river signals for each river. Default is None.

    Returns:
    --------
    xarray.Dataset

    Notes:
    ------
    This function configures the river file for a ROMS simulation. It slices the river
    data based on the specified time range, assigns temperature and salinity values,
    and optionally sets river directions and signals.

    Example:
    --------
    roms_river(ds, ds_list, temp_list, salt_list, 
               np.datetime64('2021-01-01'), np.datetime64('2021-12-31'), '1D',
               dir_list=[0.0, 45.0, -90.0], sign_list=[-1.0, 1.0, -1.0])
    """

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


def fill_with_climatology(dataarray):
    """
    Fill NaN values in a DataArray with the daily climatology average.

    Parameters:
    -----------
    dataarray : xarray.DataArray
        Input data array containing NaN values to be filled.

    Returns:
    --------
    xarray.DataArray
        Data array with NaN values replaced by the daily climatology average.

    Notes:
    ------
    This function takes a DataArray containing missing values (NaN) and fills them
    with the corresponding daily climatology average. The daily climatology average
    is computed by averaging data values for each day of the year across multiple
    years. The resulting DataArray has NaN values replaced with the computed
    climatology average.

    Example:
    --------
    climatology_filled = fill_with_climatology(data_array)
    """
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



def extend_dataset_in_time(ds, tstart, tfinal, freq='D'):
    # extend a dataset by appending another dataset to its end

    end_date = pd.to_datetime(tfinal)
    new_time_range = pd.date_range(start=tstart, end=end_date, freq=freq)

    # dataset that will be appended
    appended_ds = xr.Dataset(
            {var_name: (['time'], np.full(len(new_time_range), np.nan, dtype=ds[var_name].dtype))
            for var_name in ds.data_vars},
            coords={'time': new_time_range})
    
    extended_ds = xr.concat([ds, appended_ds], dim='time')
    return extended_ds

# dictionary describing the number of points, temperature and salinity
# that will be used to represent the rivers in the river file
# the river may be separated in 'n' fractions to avoid instabilities
# in roms, specially when the flow is large
# temp is a constant temperature set throught the run
# salt is a constant salinity set throught the run
number_of_river_points = {
                    'laplata'       : {'n':3, 'temp':17., 'salt':0},
                    'lagoapatos'    : {'n':2, 'temp':17., 'salt':0},
                    'tubarao'       : {'n':1, 'temp':20., 'salt':0},
                    'itajai'        : {'n':1, 'temp':20., 'salt':0},
                    'itapocu'       : {'n':1, 'temp':20., 'salt':0},
                    'iguape'        : {'n':1, 'temp':20., 'salt':0},
                    'mambucaba'     : {'n':1, 'temp':20., 'salt':0},
                    'macacu'        : {'n':1, 'temp':20., 'salt':0},
                    'paraibasul'    : {'n':1, 'temp':20., 'salt':0}}

# read data
roms_river_template2  = '/home/informatica/Documents/2023/roms_controller/pyroms_tools/scripts_aux/rivers_discharge/river_template.nc'
grid = xr.open_dataset('/home/informatica/model/projects/ceresIV_202309_2.012/input/roms01c_grd.nc')
z = grid.roms.s2z('rho')

roms2 = xr.open_dataset(roms_river_template2)

# getting grid positions of rivers to change vertical shape
positions =  list(zip(roms2.river_Eposition.values,
                        roms2.river_Xposition.values))
roms2.river_Vshape.load()
dicts ={i: positions[i] for i in range(len(list(zip(roms2.river_Eposition.values,
                                                    roms2.river_Xposition.values))))}

ds_river_series = grdc_prepare() # getting rivers time series

# preparing limits of dataset (extend if necessary, only forward extension in time was configured)
tstart = ds_river_series['time'].values[-1] + np.timedelta64(1,'D')  # initial time of extension
tfinal = "2023-01-01"                                                # final time of extension

# create a longer time series, extended part of data arrays are filled with nans
extended_ds = extend_dataset_in_time(ds_river_series,tstart, tfinal) 
extended_ds_cut = extended_ds.sel(time=slice(roms2.river_time.min(), roms2.river_time.max()))


# create extended and sliced ds 
for d,varname in enumerate(number_of_river_points):
    aux = fill_with_climatology(extended_ds[varname])  # just filling gaps with daily climat.
    extended_ds_cut[varname].values = aux.sel(time=slice(extended_ds_cut.time.min(),
                                                         extended_ds_cut.time.max()))
                                                        
    # creating a profile for the vertica shape of velocities in the rivers
    # the flux is stronger close to the surface
    x = z[dicts[d]].values  # stablishes the vertical domain
    x = x-x.min()

    y = np.exp(x**0.7)      # this is the shape as a function of x
    y = y/y.sum()           # guarantees that the integral is 1 (ensure mass conservation)
    roms2.river_Vshape.values[:,d]=y


# setting up time arguments in roms_river method
tstart = roms2.river_time.min()
tend   = roms2.river_time.max()
timefreq = '1D'

# preparing empty lists used in roms_river method
# river discharge fractions are used to avoid instability issues in ROMS
dslist = []  # this will have the  values of all river discharge fractions
temp   = []  # list of temperatures of each river fraction (constant throughout the time)
salt   = []  # list of salinities of each river fraction (constant throughout the time)


# loop over rivers to create lists corresponding to each river fraction
# the lists are ingested by roms_river method
for i in number_of_river_points:
    for j in range(number_of_river_points[i]['n']):
        dslist.append(extended_ds_cut[i]/number_of_river_points[i]['n'])
        temp.append(number_of_river_points[i]['temp'])
        salt.append(number_of_river_points[i]['salt'])

# create dataset 
ds_river_roms = roms_river(roms2,
                    dslist,
                    np.array(temp),
                    salt,
                    tstart,tend,timefreq,)

ds_river_roms.to_netcdf('rivers.nc')
