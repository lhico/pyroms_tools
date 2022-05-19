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




def remap_bdry(src_file, src_grd, dst_grd, src_varname, starttime,
               dst_file='', weight_file='', dmax=0,
               cdepth=0, kk=0, strtime='%Y-%m-%dT%H:%M:%S'):

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

    xrange = src_grd.xrange
    yrange = src_grd.yrange

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
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    # get missing values
    spval = src_var._FillValue
    src_var = np.squeeze(src_var)
    src_var[src_var == spval] = 1e37
    spval = 1e37
    ndim = src_var.ndim

    if ndim == 3:
        src_var = src_var[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif ndim == 2:
        src_var = src_var[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    print("dimensions:", src_var.shape, ndim)

    Cpos = 'rho'
    z = src_grd.z_t
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    local_dict = dstvar_dict[f'{src_varname}']
    wts_file = weight_file

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


    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1, 0, 0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel,
                                               nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid,
                                         dst_zcoord)

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

    # remapping
    print('remapping', dst_varname_north, 'from', src_grd.name, 'to', dst_grd.name)
    print('time =', time)

    if ndim == 3:
        # flood the grid
        print('flood the grid')
        src_varz = pyroms_toolbox.CGrid_GLORYS.flood(src_var, src_grd,
                                Cpos='t', spval=spval,
                                dmax=dmax, cdepth=cdepth, kk=kk)
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('vertical interpolation from standard z level to sigma')
        dst_var_north = z2roms(dst_varz[::-1, Mp-1:Mp, 0:Lp],
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval,
                          flood=False, irange=(0, Lp), jrange=(Mp-1 ,Mp))
        dst_var_south = z2roms(dst_varz[::-1, 0:1, :],
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval,
                          flood=False, irange=(0, Lp), jrange=(0,1))
        dst_var_east = z2roms(dst_varz[::-1, :, Lp-1:Lp],
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval,
                          flood=False, irange=(Lp-1, Lp), jrange=(0, Mp))
        dst_var_west = z2roms(dst_varz[::-1, :, 0:1],
                          dst_grdz, dst_grd, Cpos=Cpos, spval=spval,
                          flood=False, irange=(0, 1), jrange=(0, Mp))
    else:
        dst_var_north = dst_varz[-1, :]
        dst_var_south = dst_varz[0, :]
        dst_var_east = dst_varz[:, -1]
        dst_var_west = dst_varz[:, 0]

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = timeNum
    nc.variables[dst_varname_north][0] = np.squeeze(dst_var_north)
    nc.variables[dst_varname_south][0] = np.squeeze(dst_var_south)
    nc.variables[dst_varname_east][0] = np.squeeze(dst_var_east)
    nc.variables[dst_varname_west][0] = np.squeeze(dst_var_west)

    # close file
    nc.close()
    cdf.close()


def remap_bdry_uv(src_file, src_grd, dst_grd, starttime,
                  dst_fileu='', dst_filev='',
                  srcu='uo', srcv='vo',
                  strtime='%Y-%m-%dT%H:%M:%S',
                  dmax=0, cdepth=0, kk=0, wts_file=[]):

    xrange = src_grd.xrange
    yrange = src_grd.yrange


    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime, strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref
    # timeNum = timeNum + 2.5  # 5-day average

    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape

    # load var
    cdf = netCDF.Dataset(src_file)
    src_varu = cdf.variables[srcu]
    src_varv = cdf.variables[srcv]

    # -- get missing values --#
    spval = src_varu._FillValue
    src_varu = np.squeeze(src_varu)
    src_varu[src_varu == spval] = 1e37
    spval = 1e37
    ndim = src_varu.ndim

    spval = src_varv._FillValue
    src_varv = np.squeeze(src_varv)
    src_varv[src_varv == spval] = 1e37
    spval = 1e37
    ndim = src_varv.ndim

    # -- cutting data from src netdf file -- #
    src_varu = src_varu[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    src_varv = src_varv[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    print("dimensions:", src_varu.shape, ndim)

    # -- get weights file -- #
    wts_file_a = wts_file[0]
    wts_file_u = wts_file[1]
    wts_file_v = wts_file[2]

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1, 0, 0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid,
                                     dst_zcoord)

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

    # remaping
    print('remapping and rotating u and v from', src_grd.name,
                      'to', dst_grd.name)
    print('time =', time)

    # flood the grid
    print('flood the grid')
    src_uz = pyroms_toolbox.CGrid_GLORYS.flood(src_varu, src_grd, Cpos='u',
                spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
    src_vz = pyroms_toolbox.CGrid_GLORYS.flood(src_varv, src_grd, Cpos='v',
                spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_uz = pyroms.remapping.remap(src_uz, wts_file_u,
                                      spval=spval)
    dst_vz = pyroms.remapping.remap(src_vz, wts_file_v,
                                      spval=spval)

    # vertical interpolation from standard z level to sigma
    print('vertical interpolation from standard z level to sigma')
    dst_u_north = z2roms(dst_uz[::-1, Mp-2:Mp, 0:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_u_south = z2roms(dst_uz[::-1, 0:2, 0:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_u_east = z2roms(dst_uz[::-1, 0:Mp, Lp-2:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
    dst_u_west = z2roms(dst_uz[::-1, 0:Mp, 0:2],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,2), jrange=(0,Mp))

    dst_v_north = z2roms(dst_vz[::-1, Mp-2:Mp, 0:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_v_south = z2roms(dst_vz[::-1, 0:2, 0:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_v_east = z2roms(dst_vz[::-1, 0:Mp, Lp-2:Lp],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
    dst_v_west = z2roms(dst_vz[::-1, 0:Mp, 0:2],
                      dst_grdz, dst_grd, Cpos='rho', spval=spval,
                      flood=False, irange=(0,2), jrange=(0,Mp))


    # rotate u,v fields
    src_angle = src_grd.angle
    src_angle = pyroms.remapping.remap(src_angle, wts_file_a)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U_north = dst_u_north + dst_v_north*1j
    eitheta_north = np.exp(-1j*angle[:,Mp-2:Mp, 0:Lp])
    U_north = U_north * eitheta_north
    dst_u_north = np.real(U_north)
    dst_v_north = np.imag(U_north)

    U_south = dst_u_south + dst_v_south*1j
    eitheta_south = np.exp(-1j*angle[:,0:2, 0:Lp])
    U_south = U_south * eitheta_south
    dst_u_south = np.real(U_south)
    dst_v_south = np.imag(U_south)

    U_east = dst_u_east + dst_v_east*1j
    eitheta_east = np.exp(-1j*angle[:,0:Mp, Lp-2:Lp])
    U_east = U_east * eitheta_east
    dst_u_east = np.real(U_east)
    dst_v_east = np.imag(U_east)

    U_west = dst_u_west + dst_v_west*1j
    eitheta_west = np.exp(-1j*angle[:,0:Mp, 0:2])
    U_west = U_west * eitheta_west
    dst_u_west = np.real(U_west)
    dst_v_west = np.imag(U_west)

    # move back to u,v points
    dst_u_north = 0.5 * np.squeeze(dst_u_north[:,-1,:-1] + dst_u_north[:,-1,1:])
    dst_v_north = 0.5 * np.squeeze(dst_v_north[:,:-1,:] + dst_v_north[:,1:,:])
    dst_u_south = 0.5 * np.squeeze(dst_u_south[:,0,:-1] + dst_u_south[:,0,1:])
    dst_v_south = 0.5 * np.squeeze(dst_v_south[:,:-1,:] + dst_v_south[:,1:,:])
    dst_u_east = 0.5 * np.squeeze(dst_u_east[:,:,:-1] + dst_u_east[:,:,1:])
    dst_v_east = 0.5 * np.squeeze(dst_v_east[:,:-1,-1] + dst_v_east[:,1:,-1])
    dst_u_west = 0.5 * np.squeeze(dst_u_west[:,:,:-1] + dst_u_west[:,:,1:])
    dst_v_west = 0.5 * np.squeeze(dst_v_west[:,:-1,0] + dst_v_west[:,1:,0])

    if dst_u_north.ndim == 1:
        dst_u_north = dst_u_north.reshape(1,dst_u_north.size)
        dst_u_south = dst_u_south.reshape(1,dst_u_south.size)
        dst_u_east = dst_u_east.reshape(1,dst_u_east.size)
        dst_u_west = dst_u_west.reshape(1,dst_u_west.size)
        dst_v_north = dst_v_north.reshape(1,dst_v_north.size)
        dst_v_south = dst_v_south.reshape(1,dst_v_south.size)
        dst_v_east = dst_v_east.reshape(1,dst_v_east.size)
        dst_v_west = dst_v_west.reshape(1,dst_v_west.size)


    # spval
    idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
    idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
    idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
    idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
    idxu_east = np.where(dst_grd.hgrid.mask_u[:,-1] == 0)
    idxv_east = np.where(dst_grd.hgrid.mask_v[:,-1] == 0)
    idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
    idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u_north[n, idxu_north[0]] = spval
        dst_v_north[n, idxv_north[0]] = spval
        dst_u_south[n, idxu_south[0]] = spval
        dst_v_south[n, idxv_south[0]] = spval
        dst_u_east[n, idxu_east[0]] = spval
        dst_v_east[n, idxv_east[0]] = spval
        dst_u_west[n, idxu_west[0]] = spval
        dst_v_west[n, idxv_west[0]] = spval

    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] + dst_grd.vgrid.z_w[0,:,-1,1:])
    z_v_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] + dst_grd.vgrid.z_w[0,:,-2,:])
    z_u_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] + dst_grd.vgrid.z_w[0,:,0,1:])
    z_v_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] + dst_grd.vgrid.z_w[0,:,1,:])
    z_u_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:,-1] + dst_grd.vgrid.z_w[0,:,:,-2])
    z_v_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,-1] + dst_grd.vgrid.z_w[0,:,1:,-1])
    z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])
    z_v_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] + dst_grd.vgrid.z_w[0,:,1:,0])

    dst_ubar_north = np.zeros(dst_u_north.shape[1])
    dst_ubar_south = np.zeros(dst_u_south.shape[1])
    dst_ubar_east = np.zeros(dst_u_east.shape[1])
    dst_ubar_west = np.zeros(dst_u_west.shape[1])
    dst_vbar_north = np.zeros(dst_v_north.shape[1])
    dst_vbar_south = np.zeros(dst_v_south.shape[1])
    dst_vbar_east = np.zeros(dst_v_east.shape[1])
    dst_vbar_west = np.zeros(dst_v_west.shape[1])

    for i in range(dst_u_north.shape[1]):
        dst_ubar_north[i] = (dst_u_north[:,i] * np.diff(z_u_north[:,i])).sum() / -z_u_north[0,i]
        dst_ubar_south[i] = (dst_u_south[:,i] * np.diff(z_u_south[:,i])).sum() / -z_u_south[0,i]
    for i in range(dst_v_north.shape[1]):
        dst_vbar_north[i] = (dst_v_north[:,i] * np.diff(z_v_north[:,i])).sum() / -z_v_north[0,i]
        dst_vbar_south[i] = (dst_v_south[:,i] * np.diff(z_v_south[:,i])).sum() / -z_v_south[0,i]
    for j in range(dst_u_east.shape[1]):
        dst_ubar_east[j] = (dst_u_east[:,j] * np.diff(z_u_east[:,j])).sum() / -z_u_east[0,j]
        dst_ubar_west[j] = (dst_u_west[:,j] * np.diff(z_u_west[:,j])).sum() / -z_u_west[0,j]
    for j in range(dst_v_east.shape[1]):
        dst_vbar_east[j] = (dst_v_east[:,j] * np.diff(z_v_east[:,j])).sum() / -z_v_east[0,j]
        dst_vbar_west[j] = (dst_v_west[:,j] * np.diff(z_v_west[:,j])).sum() / -z_v_west[0,j]

    #mask
    dst_ubar_north = np.ma.masked_where(dst_grd.hgrid.mask_u[-1,:] == 0, dst_ubar_north)
    dst_ubar_south = np.ma.masked_where(dst_grd.hgrid.mask_u[0,:] == 0, dst_ubar_south)
    dst_ubar_east = np.ma.masked_where(dst_grd.hgrid.mask_u[:,-1] == 0, dst_ubar_east)
    dst_ubar_west = np.ma.masked_where(dst_grd.hgrid.mask_u[:,0] == 0, dst_ubar_west)
    dst_vbar_north = np.ma.masked_where(dst_grd.hgrid.mask_v[-1,:] == 0, dst_vbar_north)
    dst_vbar_south = np.ma.masked_where(dst_grd.hgrid.mask_v[0,:] == 0, dst_vbar_south)
    dst_vbar_east = np.ma.masked_where(dst_grd.hgrid.mask_v[:,-1] == 0, dst_vbar_east)
    dst_vbar_west = np.ma.masked_where(dst_grd.hgrid.mask_v[:,0] == 0, dst_vbar_west)

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = timeNum
    ncu.variables['u_north'][0] = dst_u_north
    ncu.variables['u_south'][0] = dst_u_south
    ncu.variables['u_east'][0] = dst_u_east
    ncu.variables['u_west'][0] = dst_u_west
    ncu.variables['ubar_north'][0] = dst_ubar_north
    ncu.variables['ubar_south'][0] = dst_ubar_south
    ncu.variables['ubar_east'][0] = dst_ubar_east
    ncu.variables['ubar_west'][0] = dst_ubar_west

    ncv.variables['ocean_time'][0] = timeNum
    ncv.variables['v_north'][0] = dst_v_north
    ncv.variables['v_south'][0] = dst_v_south
    ncv.variables['v_east'][0] = dst_v_east
    ncv.variables['v_west'][0] = dst_v_west
    ncv.variables['vbar_north'][0] = dst_vbar_north
    ncv.variables['vbar_south'][0] = dst_vbar_south
    ncv.variables['vbar_east'][0] = dst_vbar_east
    ncv.variables['vbar_west'][0] = dst_vbar_west

#    print dst_u.shape
#    print dst_ubar.shape
#    print dst_v.shape
#    print dst_vbar.shape

    # close file
    ncu.close()
    ncv.close()
    cdf.close()


def z2roms(varz, grdz, grd, Cpos='rho', irange=None, jrange=None, \
           spval=1e37, flood=True, dmax=0, cdepth=0, kk=0, \
           mode='linear'):
    """
    var = z2roms(var, grdz, grd)

    optional switch:
      - Cpos='rho', 'u' or 'v'       specify the C-grid position where
                                     the variable rely
      - irange                       specify grid sub-sample for i direction
      - jrange                       specify grid sub-sample for j direction
      - spval=1e37                   define spval value
      - dmax=0                       if dmax>0, maximum horizontal
                                     flooding distance
      - cdepth=0                     critical depth for flooding
                                     if depth<cdepth => no flooding
      - kk
      - mode='linear' or 'spline'    specify the type of interpolation

    Interpolate the variable from z vertical grid grdz to ROMS grid grd
    """

    varz = varz.copy()

    assert len(varz.shape) == 3, 'var must be 3D'

    if mode=='linear':
        imode=0
    elif mode=='spline':
        imode=1
    else:
        raise Warning('%s not supported, defaulting to linear' % mode)

    if Cpos == 'rho':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_r[0,:]
        mask = grd.hgrid.mask_rho
    elif Cpos == 'u':
        z = 0.5 * (grdz.vgrid.z[:,:,:-1] + grdz.vgrid.z[:,:,1:])
        depth = 0.5 * (grd.vgrid.z_r[0,:,:,:-1] + grd.vgrid.z_r[0,:,:,1:])
        mask = grd.hgrid.mask_u
    elif Cpos == 'v':
        z = 0.5 * (grdz.vgrid.z[:,:-1,:] + grdz.vgrid.z[:,1:,:])
        depth = 0.5 * (grd.vgrid.z_r[0,:,:-1,:] + grd.vgrid.z_r[0,:,1:,:])
        mask = grd.hgrid.mask_v
    elif Cpos == 'w':
        z = grdz.vgrid.z[:]
        depth = grd.vgrid.z_w[0,:]
        mask = grd.hgrid.mask_rho
    else:
        raise Warning('%s bad position. Use depth at Arakawa-C \
                             rho points instead.' % Cpos)

    nlev, Mm, Lm = varz.shape

    if depth.ndim==2:
        shape = depth.shape
        depth = depth.reshape([1,shape[0],shape[1]])

    Nm = depth.shape[0]

    if irange is None:
        irange = (0,Lm)
    else:
        assert varz.shape[2] == irange[1]-irange[0], \
               'var shape and irange must agree'

    if jrange is None:
        jrange = (0,Mm)
    else:
        assert varz.shape[1] == jrange[1]-jrange[0], \
               'var shape and jrange must agree'

    # flood varz if requested
    if flood:
        varz = pyroms.remapping.flood(varz, grdz, Cpos=Cpos, \
                 irange=irange, jrange=jrange, spval=spval, \
                 dmax=dmax, cdepth=cdepth, kk=kk)

    varz = np.concatenate((varz[0:1,:,:], varz, varz[-1:,:,:]), 0)
    z = np.concatenate((-9999*np.ones((1,z.shape[1], z.shape[2])), \
           z, \
           100*np.ones((1,z.shape[1], z.shape[2]))), 0)

    var = np.ma.zeros((Nm, Mm, Lm))
    

    for k in range(Nm):
        var[k,:,:] = pyroms._interp.xhslice(varz, \
                      z[:,jrange[0]:jrange[1], irange[0]:irange[1]], \
                      depth[k,jrange[0]:jrange[1], irange[0]:irange[1]], \
                      mask[jrange[0]:jrange[1], irange[0]:irange[1]], \
                      imode, spval)
        #mask
        var = np.ma.masked_values(var, spval, rtol=1e-5)
        #var[k,:,:] = np.ma.masked_where(mask == 0, var[k,:,:])

    return var