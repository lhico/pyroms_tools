from distutils.command.sdist import sdist
import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
from pyroms import _remapping

class nctime(object):
    pass


def make_bdry_rho_file(dst_grd, src_varname, starttime,
               dst_file='', weight_file='', strtime='%Y-%m-%dT%H:%M:%S'):

    dstvar_dict = {
      'zos': {
        'varname': 'zeta',
        'netcdf': {
            'dst_varname_north': 'zeta_north',
            'dimensions_north': ('ocean_time', 'xi_rho'),
            'long_name_north': 'free-surface north boundary condition',
            'field_north': 'zeta_north, scalar, series',
            'dst_varname_south': 'zeta_south',
            'dimensions_south': ('ocean_time', 'xi_rho'),
            'long_name_south': 'free-surface south boundary condition',
            'field_south': 'zeta_south, scalar, series',
            'dst_varname_east': 'zeta_east',
            'dimensions_east': ('ocean_time', 'eta_rho'),
            'long_name_east': 'free-surface east boundary condition',
            'field_east': 'zeta_east, scalar, series',
            'dst_varname_west': 'zeta_west',
            'dimensions_west': ('ocean_time', 'eta_rho'),
            'long_name_west': 'free-surface west boundary condition',
            'field_west': 'zeta_west, scalar, series',
            'units': 'meter',
            }
            },
      'thetao': {
        'varname': 'temp',
        'netcdf': {
            'dst_varname_north': 'temp_north',
            'dimensions_north': ('ocean_time', 's_rho', 'xi_rho'),
            'long_name_north': 'potential temperature north boundary condition',
            'field_north': 'temp_north, scalar, series',
            'dst_varname_south': 'temp_south',
            'dimensions_south': ('ocean_time', 's_rho', 'xi_rho'),
            'long_name_south': 'potential temperature south boundary condition',
            'field_south': 'temp_south, scalar, series',
            'dst_varname_east': 'temp_east',
            'dimensions_east': ('ocean_time', 's_rho', 'eta_rho'),
            'long_name_east': 'potential temperature east boundary condition',
            'field_east': 'temp_east, scalar, series',
            'dst_varname_west': 'temp_west',
            'dimensions_west': ('ocean_time', 's_rho', 'eta_rho'),
            'long_name_west': 'potential temperature west boundary condition',
            'field_west': 'temp_west, scalar, series',
            'units': 'Celsius',
        }
      },
      'so': {
        'varname': 'salt',
        'netcdf': {
            'dst_varname_north': 'salt_north',
            'dimensions_north': ('ocean_time', 's_rho', 'xi_rho'),
            'long_name_north': 'salinity north boundary condition',
            'field_north': 'salt_north, scalar, series',
            'dst_varname_south': 'salt_south',
            'dimensions_south': ('ocean_time', 's_rho', 'xi_rho'),
            'long_name_south': 'salinity south boundary condition',
            'field_south': 'salt_south, scalar, series',
            'dst_varname_east': 'salt_east',
            'dimensions_east': ('ocean_time', 's_rho', 'eta_rho'),
            'long_name_east': 'salinity east boundary condition',
            'field_east': 'salt_east, scalar, series',
            'dst_varname_west': 'salt_west',
            'dimensions_west': ('ocean_time', 's_rho', 'eta_rho'),
            'long_name_west': 'salinity west boundary condition',
            'field_west': 'salt_west, scalar, series',
            'units': 'PSU'
          }
        }
      }


    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime, strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref
    # timeNum = timeNum + 2.5  # 5-day average


    # -- read source netcdf -- #
    # cdf = netCDF.Dataset(src_file)
    # src_var = cdf.variables[src_varname]

    # get missing values
    # spval = src_var._FillValue
    # # src_var = np.squeeze(src_var)
    # # src_var[src_var == spval] = 1e37
    spval = 1e37
    # ndim = src_var.ndim

    Cpos = 'rho'
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    local_dict = dstvar_dict[f'{src_varname}']
   

    dst_varname_north = local_dict['netcdf']['dst_varname_north']
    dimensions_north = local_dict['netcdf']['dimensions_north']
    long_name_north = local_dict['netcdf']['long_name_north']
    field_north = local_dict['netcdf']['field_north']

    dst_varname_south = local_dict['netcdf']['dst_varname_south']
    dimensions_south = local_dict['netcdf']['dimensions_south']
    long_name_south = local_dict['netcdf']['long_name_south']
    field_south = local_dict['netcdf']['field_south']

    dst_varname_west = local_dict['netcdf']['dst_varname_west']
    dimensions_west = local_dict['netcdf']['dimensions_west']
    long_name_west = local_dict['netcdf']['long_name_west']
    field_west = local_dict['netcdf']['field_west']

    dst_varname_east = local_dict['netcdf']['dst_varname_east']
    dimensions_east = local_dict['netcdf']['dimensions_east']
    long_name_east = local_dict['netcdf']['long_name_east']
    field_east = local_dict['netcdf']['field_east']
    units = local_dict['netcdf']['units']


    # -- create netcdf boundary file -- #
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open boundary file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    # create variable in boudary file
    print('Creating variable', dst_varname_north)
    nc.createVariable(dst_varname_north, 'f8', dimensions_north,
                      fill_value=spval)
    nc.variables[dst_varname_north].long_name = long_name_north
    nc.variables[dst_varname_north].units = units
    nc.variables[dst_varname_north].field = field_north

    print('Creating variable', dst_varname_south)
    nc.createVariable(dst_varname_south, 'f8', dimensions_south,
                      fill_value=spval)
    nc.variables[dst_varname_south].long_name = long_name_south
    nc.variables[dst_varname_south].units = units
    nc.variables[dst_varname_south].field = field_south

    print('Creating variable', dst_varname_east)
    nc.createVariable(dst_varname_east, 'f8', dimensions_east,
                      fill_value=spval)
    nc.variables[dst_varname_east].long_name = long_name_east
    nc.variables[dst_varname_east].units = units
    nc.variables[dst_varname_east].field = field_east

    print('Creating variable', dst_varname_west)
    nc.createVariable(dst_varname_west, 'f8', dimensions_west,
                      fill_value=spval)
    nc.variables[dst_varname_west].long_name = long_name_west
    nc.variables[dst_varname_west].units = units
    nc.variables[dst_varname_west].field = field_west



    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = timeNum
    nc.variables[dst_varname_north][0] = np.nan
    nc.variables[dst_varname_south][0] = np.nan
    nc.variables[dst_varname_east][0] = np.nan
    nc.variables[dst_varname_west][0] = np.nan

    # close file
    nc.close()



def make_bdry_uv_file(dst_grd, starttime,
                  dst_fileu='', dst_filev='',
                  strtime='%Y-%m-%dT%H:%M:%S',
                  ):



    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime, strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref

    spval = 1e37


    # -- create a roms file -- #
    pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
    pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)

    # -- open destination file -- #
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    # create variable in destination file
    print('Creating variable u_north')
    ncu.createVariable('u_north', 'f8', ('ocean_time', 's_rho', 'xi_u'),
                       fill_value=spval)
    ncu.variables['u_north'].long_name = '3D u-momentum north boundary condition'
    ncu.variables['u_north'].units = 'meter second-1'
    ncu.variables['u_north'].field = 'u_north, scalar, series'
    print('Creating variable u_south')
    ncu.createVariable('u_south', 'f8', ('ocean_time', 's_rho', 'xi_u'),
                       fill_value=spval)
    ncu.variables['u_south'].long_name = '3D u-momentum south boundary condition'
    ncu.variables['u_south'].units = 'meter second-1'
    ncu.variables['u_south'].field = 'u_south, scalar, series'
    print('Creating variable u_east')
    ncu.createVariable('u_east', 'f8', ('ocean_time', 's_rho', 'eta_u'),
                       fill_value=spval)
    ncu.variables['u_east'].long_name = '3D u-momentum east boundary condition'
    ncu.variables['u_east'].units = 'meter second-1'
    ncu.variables['u_east'].field = 'u_east, scalar, series'
    print('Creating variable u_west')
    ncu.createVariable('u_west', 'f8', ('ocean_time', 's_rho', 'eta_u'),
                       fill_value=spval)
    ncu.variables['u_west'].long_name = '3D u-momentum west boundary condition'
    ncu.variables['u_west'].units = 'meter second-1'
    ncu.variables['u_west'].field = 'u_east, scalar, series'

    # create variable in destination file
    print('Creating variable ubar_north')
    ncu.createVariable('ubar_north', 'f8', ('ocean_time', 'xi_u'),
                       fill_value=spval)
    ncu.variables['ubar_north'].long_name = '2D u-momentum north boundary condition'
    ncu.variables['ubar_north'].units = 'meter second-1'
    ncu.variables['ubar_north'].field = 'ubar_north, scalar, series'
    print('Creating variable ubar_south')
    ncu.createVariable('ubar_south', 'f8', ('ocean_time', 'xi_u'),
                       fill_value=spval)
    ncu.variables['ubar_south'].long_name = '2D u-momentum south boundary condition'
    ncu.variables['ubar_south'].units = 'meter second-1'
    ncu.variables['ubar_south'].field = 'ubar_south, scalar, series'
    print('Creating variable ubar_east')
    ncu.createVariable('ubar_east', 'f8', ('ocean_time', 'eta_u'),
                       fill_value=spval)
    ncu.variables['ubar_east'].long_name = '2D u-momentum east boundary condition'
    ncu.variables['ubar_east'].units = 'meter second-1'
    ncu.variables['ubar_east'].field = 'ubar_east, scalar, series'
    print('Creating variable ubar_west')
    ncu.createVariable('ubar_west', 'f8', ('ocean_time', 'eta_u'),
                       fill_value=spval)
    ncu.variables['ubar_west'].long_name = '2D u-momentum west boundary condition'
    ncu.variables['ubar_west'].units = 'meter second-1'
    ncu.variables['ubar_west'].field = 'ubar_east, scalar, series'

    print('Creating variable v_north')
    ncv.createVariable('v_north', 'f8', ('ocean_time', 's_rho', 'xi_v'),
                       fill_value=spval)
    ncv.variables['v_north'].long_name = '3D v-momentum north boundary condition'
    ncv.variables['v_north'].units = 'meter second-1'
    ncv.variables['v_north'].field = 'v_north, scalar, series'
    print('Creating variable v_south')
    ncv.createVariable('v_south', 'f8', ('ocean_time', 's_rho', 'xi_v'),
                       fill_value=spval)
    ncv.variables['v_south'].long_name = '3D v-momentum south boundary condition'
    ncv.variables['v_south'].units = 'meter second-1'
    ncv.variables['v_south'].field = 'v_south, scalar, series'
    print('Creating variable v_east')
    ncv.createVariable('v_east', 'f8', ('ocean_time', 's_rho', 'eta_v'),
                       fill_value=spval)
    ncv.variables['v_east'].long_name = '3D v-momentum east boundary condition'
    ncv.variables['v_east'].units = 'meter second-1'
    ncv.variables['v_east'].field = 'v_east, scalar, series'
    print('Creating variable v_west')
    ncv.createVariable('v_west', 'f8', ('ocean_time', 's_rho', 'eta_v'),
                       fill_value=spval)
    ncv.variables['v_west'].long_name = '3D v-momentum west boundary condition'
    ncv.variables['v_west'].units = 'meter second-1'
    ncv.variables['v_west'].field = 'v_east, scalar, series'

    print('Creating variable vbar_north')
    ncv.createVariable('vbar_north', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_north'].long_name = '2D v-momentum north boundary condition'
    ncv.variables['vbar_north'].units = 'meter second-1'
    ncv.variables['vbar_north'].field = 'vbar_north, scalar, series'
    print('Creating variable vbar_south')
    ncv.createVariable('vbar_south', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_south'].long_name = '2D v-momentum south boundary condition'
    ncv.variables['vbar_south'].units = 'meter second-1'
    ncv.variables['vbar_south'].field = 'vbar_south, scalar, series'
    print('Creating variable vbar_east')
    ncv.createVariable('vbar_east', 'f8', ('ocean_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_east'].long_name = '2D v-momentum east boundary condition'
    ncv.variables['vbar_east'].units = 'meter second-1'
    ncv.variables['vbar_east'].field = 'vbar_east, scalar, series'
    print('Creating variable vbar_west')
    ncv.createVariable('vbar_west', 'f8', ('ocean_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_west'].long_name = '2D v-momentum west boundary condition'
    ncv.variables['vbar_west'].units = 'meter second-1'
    ncv.variables['vbar_west'].field = 'vbar_east, scalar, series'


    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = timeNum
    ncu.variables['u_north'][0] = np.nan
    ncu.variables['u_south'][0] = np.nan
    ncu.variables['u_east'][0] = np.nan
    ncu.variables['u_west'][0] = np.nan
    ncu.variables['ubar_north'][0] = np.nan
    ncu.variables['ubar_south'][0] = np.nan
    ncu.variables['ubar_east'][0] = np.nan
    ncu.variables['ubar_west'][0] = np.nan

    ncv.variables['ocean_time'][0] = timeNum
    ncv.variables['v_north'][0] = np.nan
    ncv.variables['v_south'][0] = np.nan
    ncv.variables['v_east'][0] = np.nan
    ncv.variables['v_west'][0] = np.nan
    ncv.variables['vbar_north'][0] = np.nan
    ncv.variables['vbar_south'][0] = np.nan
    ncv.variables['vbar_east'][0] = np.nan
    ncv.variables['vbar_west'][0] = np.nan

#    print dst_u.shape
#    print dst_ubar.shape
#    print dst_v.shape
#    print dst_vbar.shape

    # close file
    ncu.close()
    ncv.close()