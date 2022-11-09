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
import glob
import skfmm


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
    mask[2:-2,:-2] = 0
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

def extrapolation_nearest(x,y,var, dx, maskvalue=None):
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


        # use fast marching method to avoid nearest interpolation
        # onto all nans.
        # fasting marching method 'propagates' a boundary values
        # in space. It 'spreads' the boundary. So we use this property
        # to avoid interpolating nearest neighbor values onto
        # all grid points with nan
        mask = np.isnan(varx[k])
        try:
            dd = skfmm.distance(mask.astype(float), dx=dx)
            dd[dd>2] = np.nan
        except:
            dd = np.zeros(mask.shape)*np.nan


        print(f'nearest extrapolation: {k/n*100:0.2f}%', end='\r')
        idxnan = np.where(np.isnan(varx[k]) & ~np.isnan(dd))  # point with no data
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

    reference = 'bacia_santos_nested'
    reference = 'swatl_2022'


    dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
    dicts = dicts[reference]

    # outfile = dicts['output_file']        # output file name
    rename_coords = dicts['rename_dims']  # renaming source file coordinates
    rename_vars   = dicts['rename_vars']  # renaming sourfe file variables
    varbs = dicts['varbs_rho']            # which variables will be interpolated
    invert_depth = dicts['invert']
    zdel = dicts['delete_idepths']
    # average grid spacing in degrees this is used in the fast marching method
    # within extrapolation_nearest method

    # if you boundaries are presenting null values at ocean points this is value could be the culprit
    dx   = dicts['bdry.dxdy']  

    # 1) read data
    nc_roms_grd   = xr.open_dataset(dicts['grid_dir'])   # roms grid
    nc_aux0       = xr.open_mfdataset(dicts['src_file'], concat_dim='ocean_time', combine='nested')  # auxilary file (initial condition roms file)
    nc_out0       = xr.open_mfdataset(dicts['bdry.bdry_file'], concat_dim='ocean_time', combine='nested') # target file (boundary condition roms file)
    ncaux         = nc_aux0.copy()
    outfile       = dicts['bdry.outfile']

    nc_ini_src0   = xr.open_mfdataset(dicts['bdry.src_file'], concat_dim='time', combine='nested')  # auxilary file (initial condition roms file)

    dsaux = nc_roms_grd  #.rename({'lon_rho':'lon', 'lat_rho':'lat'})  # rename variables so xesmf understand them


    for i in range(nc_ini_src0.time.size):
        nc_ini_src = nc_ini_src0.isel(time=[i])

        # check if file was already created (in case it was, the loop starts a new iteration)
        if glob.glob(outfile % (str(nc_ini_src.time.values[0])[:19])):
            print(f'{outfile % (str(nc_ini_src.time.values[0])[:19])} already saved')
            continue

        nc_out1    = nc_out0.copy()
        nc_ini_src.load()

        nc_ini_src = nc_ini_src.rename_dims(rename_coords)
        nc_ini_src = nc_ini_src.rename_vars(rename_vars)
        print(nc_ini_src.time.values)
        # nc_ini_src = nc_ini_src.assign_coords(z=nc_ini_src.z)
        
        if horizonta_homog_fields:  # se    tting homogenous horizontal fields
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
                raise IOError("time should be a coordinate")


            if nc_ini_src[varb].values.ndim == 4:
                nc_ini_src[varb].values[0] = extrapolation_nearest(nc_ini_src.longitude.values,
                                                                nc_ini_src.latitude.values,
                                                                nc_ini_src[varb].values[0,],
                                                                dx)
                var = interpolation(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb][0], dsaux, gridtype=gtype)
                ncaux[varb].values = [var]
            elif nc_ini_src[varb].values.ndim == 3:
                shp = nc_ini_src[varb].values.shape
                nc_ini_src[varb].values = extrapolation_nearest(nc_ini_src.longitude.values,
                                                                nc_ini_src.latitude.values,
                                                                nc_ini_src[varb].values,
                                                                dx)
                var = interpolation2d(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb], dsaux, gridtype=gtype)
                var = var.squeeze()
            else:
                raise IndexError('array should be either 3D or 4D')
            ncaux[varb].values = [var]

        nc_aux1 = rotate_vector_field(ncaux)

        # this part must be improved (I'm not sure if this is the best way to calculate the integrated vertical vaelus)
        nc_aux1.ubar.values[:] = ((nc_aux1.u * nc_aux1.s_rho).sum(dim='s_rho') / nc_aux1.s_rho.sum()).values
        nc_aux1.vbar.values[:] = ((nc_aux1.v * nc_aux1.s_rho).sum(dim='s_rho') / nc_aux1.s_rho.sum()).values

        # copying boundary values from aux to the boundary file
        # (notice that 'north', 'south', 'west', 'east' refere
        # actually to top, bottom, left and right directions in the matrix)
        for varb in varbs:
            if nc_aux1[varb].ndim == 4:
                nc_out1[f'{varb}_north'].values = nc_aux1[varb].values[:,:,-1,:]
                nc_out1[f'{varb}_south'].values = nc_aux1[varb].values[:,:,0,:]
                nc_out1[f'{varb}_west'].values = nc_aux1[varb].values[:,:,:,0]
                nc_out1[f'{varb}_east'].values = nc_aux1[varb].values[:,:,:,-1]
            elif nc_aux1[varb].ndim == 3:
                nc_out1[f'{varb}_north'].values = nc_aux1[varb].values[:,-1,:]
                nc_out1[f'{varb}_south'].values = nc_aux1[varb].values[:,0,:]
                nc_out1[f'{varb}_west'].values = nc_aux1[varb].values[:,:,0]
                nc_out1[f'{varb}_east'].values = nc_aux1[varb].values[:,:,-1]

        for varb in ['ubar', 'vbar']:
            nc_out1[f'{varb}_north'].values = nc_aux1[varb].values[:,-1,:]
            nc_out1[f'{varb}_south'].values = nc_aux1[varb].values[:,0,:]
            nc_out1[f'{varb}_west'].values = nc_aux1[varb].values[:,:,0]
            nc_out1[f'{varb}_east'].values = nc_aux1[varb].values[:,:,-1]
        
        nc_out1.assign_coords(ocean_time=nc_ini_src.time.values)
        nc_out1.to_netcdf(outfile % (str(nc_ini_src.time.values[0])[:19]))


