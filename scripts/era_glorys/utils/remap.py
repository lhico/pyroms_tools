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
from pyroms import _remapping

class nctime(object):
    pass

def remap(src_file, src_grd, dst_grd, src_varname, starttime, dst_file='', strtime='%Y-%m-%dT%H:%M:%S', dmax=0, cdepth=0, kk=0, weight_file=''):

    flood_dict = {
      'flood': {
        'Cpos': 't',
        },
      'vert_interp': {
        'Cpos': 'rho',
        'flood': 'false'
        }
    }

    dstvar_dict = {
      'zos': {
        'varname': 'zeta',
        'netcdf': {
        'dimensions': ('ocean_time', 'eta_rho', 'xi_rho'),
        'long_name': 'free-surface',
        'units': 'meter',
        'field': 'free-surface, scalar, series',
        },
      },
      'thetao': {
        'varname': 'temp',
        'netcdf': {
        'dimensions': ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
        'long_name': 'potencial temperature',
        'units': 'Celsius',
        'field': 'temperature, scalar, series',
        }
      },
      'so': {
        'varname': 'salt',
        'netcdf': {
        'dimensions': ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
        'long_name': 'Salintiy',
        'units': 'psu',
        'field': 'salinity, scalar, series',
        }
      }
    }


    xrange=src_grd.xrange
    yrange=src_grd.yrange
    # xrange=(200, 500)
    # yrange=(200, 500)

    # -- create time attributes -- #
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    ref = datetime(1900, 1, 1, 0, 0, 0)
    ref = date2num(ref)

    time = datetime.strptime(starttime,strtime)
    timeNum = date2num(time)
    timeNum = timeNum - ref
    timeNum = timeNum + 0.5 # 1-day average

    # -- read source netcdf -- #
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    # get missing values
    # spval = 1.e37
    spval = src_var._FillValue
    add_offset = src_var.add_offset
    scale_factor = src_var.scale_factor
    src_var = np.squeeze(src_var)
    # raise ValueError()
    src_var[src_var == spval] =  1e37
    spval= 1e37
    # spval = spval*scale_factor
    ndim = src_var.ndim


    if ndim == 3:
        src_var = src_var[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        # src_var = src_var[:,np.r_[ystart:np.size(src_var,1),-1],:]
    elif ndim == 2:
        src_var = src_var[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        # src_var = src_var[np.r_[ystart:np.size(src_var,0),-1],:]
    print("dimensions:", src_var.shape, ndim)


    # if src_varname == 'sossheig':
    Bpos = 't'
    Cpos = 'rho'
    z = src_grd.z_t
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    local_dict = dstvar_dict[f'{src_varname}']
    wts_file = weight_file
    dst_varname = local_dict['varname']
    dimensions = local_dict['netcdf']['dimensions']
    long_name = local_dict['netcdf']['long_name']
    units =  local_dict['netcdf']['units']
    field =  local_dict['netcdf']['field']



    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # -- create a roms file -- #
    print('Creating variable', dst_varname)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field

    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, \
              'to', dst_grd.name)
    print('time =', time)



    if ndim == 3:
        # flood the grid
        print('flood the grid')
        fdict = flood_dict['flood']
        src_varz = pyroms_toolbox.CGrid_GLORYS.flood(src_var, src_grd,
            cdepth=cdepth, dmax=0, kk=0, spval=spval,**fdict)
        # print('flooded the grid', src_varz[:,-1,189])
        # print('flooded the grid', src_varz[:,-1,277])
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, \
                                          spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('vertical interpolation from standard z level to sigma')
        # there is a quickfix in grids.Agrid.A2C

        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, spval=spval, flood=False,
                          Cpos=Cpos)
    else:
        dst_var = dst_varz

    # write data in destination file
    # print('write data in destination file')
    nc.variables['ocean_time'][0] = timeNum
    nc.variables[dst_varname][0] = dst_var

    # close destination file
    nc.close()


    return dst_varz


def remap_uv(src_file, src_grd, dst_grd, starttime,
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
    timeNum = timeNum + 0.5 # 1-day average

    # -- read source netcdf -- #
    cdf = netCDF.Dataset(src_file)
    src_varu = cdf.variables[srcu]
    src_varv = cdf.variables[srcv]
    print("dims", src_varu.dimensions, src_varv.dimensions)

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

    # -- build intermediate zgrid -- #
    z = src_grd.z_t
    zlevel = -z[::-1,0,0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

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

    # -- remapping -- #
    print('remapping and rotating u and v from', src_grd.name,
          'to', dst_grd.name)
    print('time =', time)

    # flood the grid
    print('flood the grid', src_varu.shape)
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
    dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz,
                        dst_grd, Cpos='rho', spval=spval, flood=False)
    dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz,
                        dst_grd, Cpos='rho', spval=spval, flood=False)

    # rotate u,v fields
    src_angle = src_grd.angle
    src_angle = pyroms.remapping.remap(src_angle, wts_file_a)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U = dst_u + dst_v*1j
    eitheta = np.exp(-1j*angle[:,:,:])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)

    # move back to u,v points
    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

    # spval
    idxu = np.where(dst_grd.hgrid.mask_u == 0)
    idxv = np.where(dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n,idxu[0], idxu[1]] = spval
        dst_v[n,idxv[0], idxv[1]] = spval

    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])

    dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
    dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

    for i in range(dst_ubar.shape[1]):
        for j in range(dst_ubar.shape[0]):
            dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

    for i in range(dst_vbar.shape[1]):
        for j in range(dst_vbar.shape[0]):
            dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

    # spval
    dst_ubar[idxu[0], idxu[1]] = spval
    dst_vbar[idxv[0], idxv[1]] = spval

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = timeNum
    ncu.variables['u'][0] = dst_u
    ncu.variables['ubar'][0] = dst_ubar

    ncv.variables['ocean_time'][0] = timeNum
    ncv.variables['v'][0] = dst_v
    ncv.variables['vbar'][0] = dst_vbar

    print(dst_u.shape)
    print(dst_ubar.shape)
    print(dst_v.shape)
    print(dst_vbar.shape)

    # close destination file
    ncu.close()
    ncv.close()
    cdf.close()
