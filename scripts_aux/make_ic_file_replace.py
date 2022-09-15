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

def interpolation(fpath, nc_roms_grd, source_grid, target_grid):

    # TODO fix up the vertical axis
    nc0 = xr.open_dataset(fpath)

    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho)
    nc_roms_grd = nc_roms_grd.assign_coords(s_w = nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values


    # define vertical coordinates
    z,_ = compute_depth_layers(nc_roms_grd) 
    z = z.transpose(*('s_rho','eta_rho','xi_rho'), transpose_coords=False)
    z = z.values

    #  horizontal interpolation
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(source_grid)  # interpolating

    interpvarb = np.zeros(z.shape)

    mask = nc_roms_grd.mask_rho.values
    ind = np.where(mask!=0)

    for j,i in zip(ind[0], ind[1]):
        print(j,i)
        f = interpolate.interp1d(-interpolated.z.values,
                                  interpolated[:,j,i].values,
                                  bounds_error=False,
                                  fill_value='extrapolate',
                                  kind='slinear')
        interpvarb[:,j,i] = f(z[::,j,i])
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





reference = 'swatl_2022'
dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
dicts = dicts[reference]

outfile = dicts['output_file']
rename_coords = dicts['rename_dims']
rename_vars   = dicts['rename_vars']
varbs = dicts['varbs_rho']
invert_depth = dicts['invert']
zdel = dicts['delete_idepths']

# 1) read data
nc_roms_grd  = xr.open_dataset(dicts['grid_dir'])  # roms grid
nc_ini_src   = xr.open_dataset(dicts['ic_file'])   # initial conditions sourec
nc_out0       = xr.open_dataset(dicts['src_file']) # target file (we will replace some of its variables)
nc_out = nc_out0.copy()

ds_out = nc_roms_grd.rename({'lon_rho':'lon', 'lat_rho':'lat'})  # rename variables so xesmf understand them

nc_ini_src = nc_ini_src.rename_dims(rename_coords)
nc_ini_src = nc_ini_src.rename_vars(rename_vars)
# nc_ini_src = nc_ini_src.assign_coords(z=nc_ini_src.z)

# attention
zsel = np.arange(nc_ini_src.z.values.size)
zsel = np.delete(zsel, zdel)
nc_ini_src = nc_ini_src.isel(z=zsel)


# 2) linear interpolation from source to roms grid (horizontal) -> roms_aux
# the extrapolation procedure is not generic!!
for varb in varbs:
    if nc_ini_src[varb].values.ndim == 4:
        nc_ini_src[varb].values[0] = extrapolation_nearest(nc_ini_src.longitude.values,
                                                           nc_ini_src.latitude.values,
                                                           nc_ini_src[varb].values[0,])
        var = interpolation(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb][0], ds_out)
    elif nc_ini_src[varb].values.ndim == 3:
        nc_ini_src[varb].values = extrapolation_nearest(nc_ini_src.longitude.values,
                                                        nc_ini_src.latitude.values,
                                                        nc_ini_src[varb].values)
        var = interpolation(dicts['grid_dir'], nc_roms_grd, nc_ini_src[varb], ds_out)
    else:
        raise IndexError('array should be either 3D or 4D')
    nc_out[varb].values = [var]

nc_out.u.values[:] = 0 
nc_out.v.values[:] = 0 
nc_out.ubar.values[:] = 0 
nc_out.vbar.values[:] = 0 

os.system(f'rm {outfile}')
nc_out.to_netcdf(outfile)




# j=184
# i=88
# z,_ = compute_depth_layers(nc_out)
# x, y = z[j,i].lon_rho.values, z[j,i].lat_rho.values


# aux = nc_ini_src.assign_coords(x=nc_ini_src.lon[0,:].values, y= nc_ini_src.lat[:,0].values )
# aux1 = aux.sel(x=x, y=y, method='nearest')

# plt.plot(nc_out['temp'][0,:,j,i], z[j,i], marker='x')
# plt.plot(aux1['temp'], -aux1.z, marker='x')