from typing import IO
import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
from datetime import datetime
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox

class nctime(object):
    pass


def create_scalar_fields_3D(dst_grd, src_varname, starttime, 
          dst_file='', strtime='%Y-%m-%dT%H:%M:%S',
          spval = np.nan):

    scalar_dict = {
      'zos': {
        'varname': 'zeta',
        'netcdf': {
        'dimensions': ('ocean_time', 'eta_rho', 'xi_rho'),
        'long_name': 'free-surface',
        'units': 'meter',
        'field': 'free-surface, scalar, series'},
        },
      'thetao': {
        'varname': 'temp',
        'netcdf': {
            'dimensions': ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
            'long_name': 'potencial temperature',
            'units': 'Celsius',
            'field': 'temperature, scalar, series',}
            },
      'so': {
        'varname': 'salt',
        'netcdf': {
            'dimensions': ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
            'long_name': 'Salintiy',
            'units': 'psu',
            'field': 'salinity, scalar, series',}
          }   
        }

    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime,strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref
    timeNum = timeNum # + 0.5 # 1-day average

    # Mp, Lp = dst_grd.hgrid.mask_rho.shape
    local_dict = scalar_dict[f'{src_varname}']
    # wts_file = weight_file
    dst_varname = local_dict['varname']
    dimensions = local_dict['netcdf']['dimensions']
    long_name = local_dict['netcdf']['long_name']
    units =  local_dict['netcdf']['units']
    field =  local_dict['netcdf']['field']


    # -- create a roms file -- #
    print('Creating variable', dst_varname)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field


    nc.variables['ocean_time'][0] = timeNum
    nc.variables[dst_varname][0] = np.nan

    # close destination file
    nc.close()



def create_vector_fields(dst_grd, starttime,
             dst_fileu='', dst_filev='',
             strtime='%Y-%m-%dT%H:%M:%S',
             spval=np.nan):

    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime, strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref
    timeNum = timeNum #+ 0.5 # 1-day average

  
    # -- create a roms file -- #
    pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
    pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)

    # -- open destination file -- #
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    # -- create variable in destination file -- #
    print('Creating variable u')
    ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'),
                       fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'

    print('Creating variable ubar')
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'),
                       fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'

    print('Creating variable v')
    ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'),
                       fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    print('Creating variable vbar')
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'),
                       fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'
  

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = timeNum
    ncu.variables['u'][0] = 0.
    ncu.variables['ubar'][0] = 0.

    ncv.variables['ocean_time'][0] = timeNum
    ncv.variables['v'][0] = 0.
    ncv.variables['vbar'][0] = 0.

    # close destination file
    ncu.close()
    ncv.close()