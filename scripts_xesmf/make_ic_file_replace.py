from typing import IO
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
from scipy import interpolate
from utils import utils as ut
import os
from scipy.spatial import cKDTree


def compute_depth_layers(ds, hmin=0.1):
    """ compute depths of ROMS vertical levels (Vtransform = 2) """
    
    # compute vertical transformation functional
    S_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    S_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
    
    # compute depth of rho (layers) and w (interfaces) points
    z_rho = ds.h * S_rho
    z_w = ds.h * S_w
    
    return z_rho, z_w

def interpolation(fpath, nc_roms_grd, source_grid, target_grid, gridtype='rho'):
    target_grid = target_grid.copy()

    coords_rename = {'roms2data':
                        {
                            'rho': {'lon_rho':'lon', 'lat_rho':'lat'},
                            'u'  : {'lon_u':'lon', 'lat_u':'lat'},
                            'v'  : {'lon_v':'lon', 'lat_v':'lat'},
                        },
                     'data2roms':
                       {
                            'rho': {'lon':'lon_rho', 'lat':'lat_rho'},
                            'u'  : {'lon':'lon_u', 'lat':'lat_u'},
                            'v'  : {'lon':'lon_v', 'lat':'lat_v'},
                       }

                    }

    # TODO fix up the vertical axis
    nc0 = xr.open_dataset(fpath)

    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho)
    nc_roms_grd = nc_roms_grd.assign_coords(s_w = nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values

    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])

    # define vertical coordinates
    z,_ = compute_depth_layers(nc_roms_grd) 
    z = z.transpose(*('s_rho','eta_rho','xi_rho'), transpose_coords=False)
    z = z.values

    if gtype == 'u':
        z = z[:,:,:-1]
    elif gtype == 'v':
        z = z[:,:-1,:]


    #  horizontal interpolation
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(source_grid)  # interpolating

    interpvarb = np.zeros(z.shape)

    mask = nc_roms_grd[f'mask_{gridtype}'].values
    ind = np.where(mask!=0)

    for j,i in zip(ind[0], ind[1]):
        print(j,i, end='\r')
        f = interpolate.interp1d(-interpolated.z.values,
                                  interpolated[:,j,i].values,
                                  bounds_error=False,
                                  fill_value='extrapolate',
                                  kind='slinear')
        interpvarb[:,j,i] = f(z[::,j,i])
    return interpvarb


def interpolation2d(fpath, nc_roms_grd, source_grid, target_grid, gridtype='rho'):
    target_grid = target_grid.copy()

    coords_rename = {'roms2data':
                        {
                            'rho': {'lon_rho':'lon', 'lat_rho':'lat'},
                            'u'  : {'lon_u':'lon', 'lat_u':'lat'},
                            'v'  : {'lon_v':'lon', 'lat_v':'lat'},
                        },
                     'data2roms':
                       {
                            'rho': {'lon':'lon_rho', 'lat':'lat_rho'},
                            'u'  : {'lon':'lon_u', 'lat':'lat_u'},
                            'v'  : {'lon':'lon_v', 'lat':'lat_v'},
                       }

                    }

    # TODO fix up the vertical axis
    nc0 = xr.open_dataset(fpath)

    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho)
    nc_roms_grd = nc_roms_grd.assign_coords(s_w = nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values

    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])

    # define vertical coordinates
    z,_ = compute_depth_layers(nc_roms_grd) 
    z = z.transpose(*('s_rho','eta_rho','xi_rho'), transpose_coords=False)
    z = z.values


    #  horizontal interpolation
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(source_grid)  # interpolating

    interpvarb = interpolated.values

    # mask = nc_roms_grd[f'mask_{gridtype}'].values
    # ind = np.where(mask!=0)

    # for j,i in zip(ind[0], ind[1]):
    #     print(j,i, end='\r')
    #     f = interpolate.interp1d(-interpolated.z.values,
    #                               interpolated[:,j,i].values,
    #                               bounds_error=False,
    #                               fill_value='extrapolate',
    #                               kind='slinear')
    #     interpvarb[:,j,i] = f(z[::,j,i])
    return interpvarb

def extrapolation_nearest(x,y,var, maskvalue=None):
    "nearest-neighbor extrapolation with cKDTree method"
    
    varx = var.copy()
    # mask values to nan (if varx does not use nan as mask already)
    if maskvalue is not None:
        varx[var==maskvalue] = np.nan
    
    varz = varx.copy()  # output variable

    # if coordinates are not 2D (we make them 2D)
    if (x.ndim == 1):
        x, y = np.meshgrid(x,y)


    n = var.shape[0]  # used to print percentage
    for k in range(var.shape[0]):
        # we will look for poinst with data and extrapolate their value
        # to nearest neighbors with no data

        print(f'nearest extrapolation: {k/n*100:0.2f}%', end='\r')
        idxnan = np.where(np.isnan(varx[k]))  # point with no data
        idx    = np.where(~np.isnan(varx[k])) # point with data

        # set up arrays to be used in CKDTree
        wet = np.zeros((len(idx[0]),2)).astype(int)
        dry = np.zeros((len(idxnan[0]),2)).astype(int)
        wet[:,0] = idx[0].astype(int)
        wet[:,1] = idx[1].astype(int)
        dry[:,0] = idxnan[0].astype(int)
        dry[:,1] = idxnan[1].astype(int)
        xwet = x[wet[:,0], wet[:,1]]
        ywet = y[wet[:,0], wet[:,1]]
        xdry = x[dry[:,0], dry[:,1]]
        ydry = y[dry[:,0], dry[:,1]]

        xy = np.array([xdry, ydry])
        tree = cKDTree(np.c_[xwet, ywet])  # prepare kdtree
        _,ii = tree.query(xy.T, k=1)       # query at nearest neighbor only
        if wet.shape[0] == 0: 
            pass
        else:
            # ii-th wet  onto dry points - shape of ii and dry must be the same
            varz[k, dry[:,0], dry[:,1]] = varx[k, wet[ii,0], wet[ii,1]]

    return varz



def rotate_vector_field(nc):
    """
    Rotates vector fields in roms ic file.
    1) regrid variables in u and v coordinates to rho coordinates
    2) rotate vectors
    3) regridded and rotated variables are  transformed back to u and v coordinates


    Args:
        nc (_type_): _description_

    Returns:
        _type_: _description_
    """    
    nc_out = nc.copy()
    nc_out.load()

    # -- regrid variables in u and v coordinates to rho coordinates -- #
    target_grid = nc_out.rename({'lon_rho': 'lon', 'lat_rho':'lat'})
    source_grid = nc_out.rename({'lon_u': 'lon', 'lat_u':'lat'})
    
    #  horizontal interpolation
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interp_u = regridder(source_grid)  # interpolating
    
    target_grid = nc_out.rename({'lon_rho': 'lon', 'lat_rho':'lat'})
    source_grid = nc_out.rename({'lon_v': 'lon', 'lat_v':'lat'})
    #  horizontal interpolation
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interp_v = regridder(source_grid)  # interpolating

    # -- rotate vectors -- #
    rot = nc_out.angle
    for var in [['u','v'],['ubar','vbar']]:
        u = interp_u[var[0]]
        v = interp_v[var[1]]

        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v,u)

        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)
        interp_u[var[0]] = u1
        interp_v[var[1]] = v1

    # -- regridded and rotated variables are  transformed back to u and v coordinates -- #
    target_grid = nc_out.rename({'lon_u': 'lon', 'lat_u':'lat'})
    source_grid = interp_u
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    regrid_u  = regridder(source_grid)

    target_grid = nc_out.rename({'lon_v': 'lon', 'lat_v':'lat'})
    source_grid = interp_v
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    regrid_v  = regridder(source_grid)

    for var in [['u','v'],['ubar','vbar']]:
        nc_out[var[0]]= regrid_u[var[0]]
        nc_out[var[1]]= regrid_v[var[1]]

    return nc_out



if __name__ == '__main__':

    # set this true when testing for horizontal homogenous fields
    # it will use the average of initial conditions source
    horizonta_homog_fields = False

    reference = 'pbs_202109_glorys'
  


    dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
    dicts = dicts[reference]

    outfile = dicts['output_file']        # output file name
    rename_coords = dicts['rename_dims']  # renaming source file coordinates
    rename_vars   = dicts['rename_vars']  # renaming sourfe file variables
    varbs = dicts['varbs_rho']            # which variables will be interpolated
    invert_depth = dicts['invert']
    zdel = dicts['delete_idepths']

    # 1) read data
    nc_roms_grd  = xr.open_dataset(dicts['grid_dir'])  # roms grid
    nc_ini_src   = xr.open_dataset(dicts['ic_file'])   # initial conditions sourec
    nc_out0       = xr.open_dataset(dicts['src_file']) # target file (we will replace some of its variables)
    nc_out = nc_out0.copy()

    ds_out = nc_roms_grd  #.rename({'lon_rho':'lon', 'lat_rho':'lat'})  # rename variables so xesmf understand them

    nc_ini_src = nc_ini_src.rename_dims(rename_coords)
    nc_ini_src = nc_ini_src.rename_vars(rename_vars)
    print(nc_ini_src.time.values)
    # nc_ini_src = nc_ini_src.assign_coords(z=nc_ini_src.z)

    
    if horizonta_homog_fields:  # setting homogenous horizontal fields
        nc = nc_ini_src.mean(dim=['lon','lat'])
        for i in varbs:
            nc_ini_src[i].values[:] = nc[i].values[:,None,None] 
        outfile = outfile[:-3] + '_hor_homog.nc'

    # attention
    zsel = np.arange(nc_ini_src.z.values.size)
    zsel = np.delete(zsel, zdel)
    nc_ini_src = nc_ini_src.isel(z=zsel)


    # 2) linear interpolation from source to roms grid (horizontal) -> roms_aux
    # salinity and temperature must have 4 dimension (time, depth, lat,)
    for varb in varbs:

        print(f'\n# -- Interpolating {varb} --#\n')
        
        gtype = 'rho'
        if varb in ['u']:
            gtype = 'u'
        elif varb in ['v']:
            gtype = 'v' 
        
        print(f'grid type: {gtype}')
        print(f'{varb} shape: {nc_ini_src[varb].values.ndim}')
        
        if "time" not in nc_ini_src[varb].coords:
            raise IOError("""time should be a coordinate. If there is a time coordinate with a different
            name, please rename it to 'time'""")


        if nc_ini_src[varb].values.ndim == 4:
            nc_ini_src[varb].values[0] = extrapolation_nearest(nc_ini_src.longitude.values,
                                                            nc_ini_src.latitude.values,
                                                            nc_ini_src[varb].values[0,])
            var = interpolation(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb][0], ds_out, gridtype=gtype)

        elif nc_ini_src[varb].values.ndim == 3:
            shp = nc_ini_src[varb].values.shape
            nc_ini_src[varb].values = extrapolation_nearest(nc_ini_src.longitude.values,
                                                            nc_ini_src.latitude.values,
                                                            nc_ini_src[varb].values)
            var = interpolation2d(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb], ds_out, gridtype=gtype)
            var = var.squeeze()
        else:
            raise IndexError('array should be either 3D or 4D')
        nc_out[varb].values = [var]

    nc_out1 = rotate_vector_field(nc_out)

    
    # nc_out.u.values[:] = 0 
    # nc_out.v.values[:] = 0 
    nc_out1.ubar.values[:] = ((nc_out1.u * nc_out1.s_rho).sum(dim='s_rho') / nc_out1.s_rho.sum()).values
    nc_out1.vbar.values[:] = ((nc_out1.v * nc_out1.s_rho).sum(dim='s_rho') / nc_out1.s_rho.sum()).values

    os.system(f'rm {outfile}')
    nc_out1.to_netcdf(outfile)

