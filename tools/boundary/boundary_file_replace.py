import sys
from typing import IO
import numpy as np
import xarray as xr
import xesmf as xe
from scipy import interpolate
from utils import utils as ut
from scipy.spatial import cKDTree
import glob
import skfmm
from netCDF4 import num2date,date2num
import pandas as pd
from datetime import datetime as dtt
import argparse
import yaml



def compute_depth_layers(ds, hmin=0.1):
    """Compute depths of ROMS vertical levels (Vtransform = 2)."""
    S_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    S_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
    
    z_rho = ds.h * S_rho
    z_w = ds.h * S_w
    
    return z_rho, z_w

def interpolation(fpath, nc_roms_grd, source_grid, target_grid, gridtype='rho'):
    target_grid = target_grid.copy()

    coords_rename = {
        'roms2data': {
            'rho': {'lon_rho': 'lon', 'lat_rho': 'lat'},
            'u': {'lon_u': 'lon', 'lat_u': 'lat'},
            'v': {'lon_v': 'lon', 'lat_v': 'lat'},
        },
        'data2roms': {
            'rho': {'lon': 'lon_rho', 'lat': 'lat_rho'},
            'u': {'lon': 'lon_u', 'lat': 'lat_u'},
            'v': {'lon': 'lon_v', 'lat': 'lat_v'},
        }
    }

    nc0 = xr.open_dataset(fpath)

    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho, s_w=nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values

    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])

    z, _ = compute_depth_layers(nc_roms_grd)
    z = z.transpose('s_rho', 'eta_rho', 'xi_rho', transpose_coords=False).values

    if gridtype == 'u':
        z = z[:, :, :-1]
    elif gridtype == 'v':
        z = z[:, :-1, :]

    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(source_grid)

    interpvarb = np.zeros(z.shape)
    mask = nc_roms_grd[f'mask_{gridtype}'].values
    mask[2:-2, :-2] = 0
    ind = np.where(mask != 0)

    for j, i in zip(ind[0], ind[1]):
        f = interpolate.interp1d(-interpolated.z.values, interpolated[:, j, i].values,
                                 bounds_error=False, fill_value='extrapolate', kind='slinear')
        interpvarb[:, j, i] = f(z[:, j, i])
    
    return interpvarb

def interpolation2d(fpath, nc_roms_grd, source_grid, target_grid, gridtype='rho'):
    target_grid = target_grid.copy()

    coords_rename = {
        'roms2data': {
            'rho': {'lon_rho': 'lon', 'lat_rho': 'lat'},
            'u': {'lon_u': 'lon', 'lat_u': 'lat'},
            'v': {'lon_v': 'lon', 'lat_v': 'lat'},
        },
        'data2roms': {
            'rho': {'lon': 'lon_rho', 'lat': 'lat_rho'},
            'u': {'lon': 'lon_u', 'lat': 'lat_u'},
            'v': {'lon': 'lon_v', 'lat': 'lat_v'},
        }
    }

    nc0 = xr.open_dataset(fpath)

    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho, s_w=nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values

    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])

    z, _ = compute_depth_layers(nc_roms_grd)
    z = z.transpose('s_rho', 'eta_rho', 'xi_rho', transpose_coords=False).values

    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(source_grid)

    return interpolated.values

def extrapolation_nearest(x, y, var, dx, maskvalue=None):
    """Nearest-neighbor extrapolation with cKDTree method."""
    varx = var.copy()
    if maskvalue is not None:
        varx[var == maskvalue] = np.nan

    varz = varx.copy()

    if x.ndim == 1:
        x, y = np.meshgrid(x, y)

    n = var.shape[0]
    for k in range(n):
        mask = np.isnan(varx[k])
        try:
            dd = skfmm.distance(mask.astype(float), dx=dx)
            dd[dd > 2] = np.nan
        except:
            dd = np.zeros(mask.shape) * np.nan

        idxnan = np.where(np.isnan(varx[k]) & ~np.isnan(dd))
        idx = np.where(~np.isnan(varx[k]))

        wet = np.column_stack((idx[0], idx[1]))
        dry = np.column_stack((idxnan[0], idxnan[1]))
        xwet = x[wet[:, 0], wet[:, 1]]
        ywet = y[wet[:, 0], wet[:, 1]]
        xdry = x[dry[:, 0], dry[:, 1]]
        ydry = y[dry[:, 0], dry[:, 1]]

        tree = cKDTree(np.c_[xwet, ywet])
        _, ii = tree.query(np.c_[xdry, ydry], k=1)
        if wet.shape[0] > 0:
            varz[k, dry[:, 0], dry[:, 1]] = varx[k, wet[ii, 0], wet[ii, 1]]

    return varz

def rotate_vector_field(nc):
    """Rotates vector fields in ROMS IC file."""
    nc_out = nc.copy()
    nc_out.load()

    target_grid = nc_out.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
    source_grid = nc_out.rename({'lon_u': 'lon', 'lat_u': 'lat'})
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interp_u = regridder(source_grid)

    source_grid = nc_out.rename({'lon_v': 'lon', 'lat_v': 'lat'})
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    interp_v = regridder(source_grid)

    rot = nc_out.angle
    for var in [['u', 'v'], ['ubar', 'vbar']]:
        u = interp_u[var[0]]
        v = interp_v[var[1]]

        mag = np.sqrt(u**2 + v**2)
        angle = np.arctan2(v, u)

        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)
        interp_u[var[0]] = u1
        interp_v[var[1]] = v1

    target_grid = nc_out.rename({'lon_u': 'lon', 'lat_u': 'lat'})
    regridder = xe.Regridder(interp_u, target_grid, 'bilinear', extrap_method='nearest_s2d')
    regrid_u = regridder(interp_u)

    target_grid = nc_out.rename({'lon_v': 'lon', 'lat_v': 'lat'})
    regridder = xe.Regridder(interp_v, target_grid, 'bilinear', extrap_method='nearest_s2d')
    regrid_v = regridder(interp_v)

    for var in [['u', 'v'], ['ubar', 'vbar']]:
        nc_out[var[0]] = regrid_u[var[0]]
        nc_out[var[1]] = regrid_v[var[1]]

    return nc_out


def generate_output_file(outfile, tref, i):
    return outfile + tref[i].strftime('%Y-%m-%dT%H:%M:%S')

def load_datasets(dicts):
    nc_roms_grd = xr.open_dataset(dicts['grid']['grid'])
    nc_aux0 = xr.open_mfdataset(dicts['bndry']['ic_file'], concat_dim='ocean_time', combine='nested')
    nc_out0 = xr.open_mfdataset(dicts['bndry']['bdry_template'], concat_dim='ocean_time', combine='nested', decode_times=False)
    nc_ini_src0 = xr.open_mfdataset(dicts['bndry']['source_file'], concat_dim='time', combine='nested', decode_times=False)
    return nc_roms_grd, nc_aux0, nc_out0, nc_ini_src0


def slice_time(nc_ini_src_, tstart, tfinal):
    tstart_num = date2num(dtt.strptime(tstart, '%Y-%m-%dT%H:%M:%S'), nc_ini_src_.time.units)
    tfinal_num = date2num(dtt.strptime(tfinal, '%Y-%m-%dT%H:%M:%S'), nc_ini_src_.time.units)
    return nc_ini_src_.sel(time=slice(tstart_num, tfinal_num + 23))

def interpolate_and_rotate(nc_ini_src, nc_aux0, dsaux, varbs, dx, dicts, nc_roms_grd):
    for varb in varbs:
        gtype = 'rho' if varb not in ['u', 'v'] else varb
        if nc_ini_src[varb].values.ndim == 4:
            nc_ini_src[varb].values[0] = extrapolation_nearest(nc_ini_src.longitude.values,
                                                               nc_ini_src.latitude.values,
                                                               nc_ini_src[varb].values[0], dx)
            var = interpolation(dicts['grid']['grid'], nc_roms_grd, nc_ini_src[varb][0], dsaux, gridtype=gtype)
            nc_aux0[varb].values = [var]
        elif nc_ini_src[varb].values.ndim == 3:
            nc_ini_src[varb].values = extrapolation_nearest(nc_ini_src.longitude.values,
                                                            nc_ini_src.latitude.values,
                                                            nc_ini_src[varb].values, dx)
            var = interpolation2d(dicts['grid']['grid'], nc_roms_grd, nc_ini_src[varb], dsaux, gridtype=gtype).squeeze()
            nc_aux0[varb].values = [var]
    nc_aux1 = rotate_vector_field(nc_aux0)
    nc_aux1.ubar.values[:] = ((nc_aux1.u * nc_aux1.s_rho).sum(dim='s_rho') / nc_aux1.s_rho.sum()).values
    nc_aux1.vbar.values[:] = ((nc_aux1.v * nc_aux1.s_rho).sum(dim='s_rho') / nc_aux1.s_rho.sum()).values
    return nc_aux1

def copy_boundary_values(nc_aux1, nc_out1, varbs):
    for varb in varbs:
        if nc_aux1[varb].ndim == 4:
            for direction in ['north', 'south', 'west', 'east']:
                nc_out1[f'{varb}_{direction}'].values = nc_aux1[varb].values[:, :, -1, :] if direction == 'north' else \
                                                        nc_aux1[varb].values[:, :, 0, :] if direction == 'south' else \
                                                        nc_aux1[varb].values[:, :, :, 0] if direction == 'west' else \
                                                        nc_aux1[varb].values[:, :, :, -1]
                nc_out1[f'{varb}_{direction}'] = nc_out1[f'{varb}_{direction}'].bfill(dim='s_rho')
        elif nc_aux1[varb].ndim == 3:
            for direction in ['north', 'south', 'west', 'east']:
                nc_out1[f'{varb}_{direction}'].values = nc_aux1[varb].values[:, -1, :] if direction == 'north' else \
                                                        nc_aux1[varb].values[:, 0, :] if direction == 'south' else \
                                                        nc_aux1[varb].values[:, :, 0] if direction == 'west' else \
                                                        nc_aux1[varb].values[:, :, -1]
    for varb in ['ubar', 'vbar']:
        for direction in ['north', 'south', 'west', 'east']:
            nc_out1[f'{varb}_{direction}'].values = nc_aux1[varb].values[:, -1, :] if direction == 'north' else \
                                                    nc_aux1[varb].values[:, 0, :] if direction == 'south' else \
                                                    nc_aux1[varb].values[:, :, 0] if direction == 'west' else \
                                                    nc_aux1[varb].values[:, :, -1]

def main():
    parser = argparse.ArgumentParser(description='Process boundary conditions.')
    parser.add_argument('--config', type=str, required=True, help='Path to the YAML configuration file')
    args = parser.parse_args()

    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)

    reference = 'ceresIV_2.012'
    dicts = config[reference]

    outfile = dicts['bndry']['output_file']
    rename_coords = dicts['bndry']['rename_dims']
    rename_vars = dicts['bndry']['map_vars']
    varbs = dicts['bndry']['varbs_rho']
    invert_depth = dicts['bndry']['invert']
    zdel = dicts['bndry']['delete_idepths']
    dx = dicts['bndry']['dxdy']
    horizontal_homog_fields = dicts['bndry']['horizontal_homog_fields']
    tstart, tfinal = dicts['bndry']['time_slice']

    nc_roms_grd, nc_aux0, nc_out0, nc_ini_src_ = load_datasets(dicts)
    nc_ini_src0 = slice_time(nc_ini_src_, tstart, tfinal)

    time0 = num2date(nc_ini_src0.time[0].values, nc_ini_src0.time.attrs['units'])
    tref = pd.date_range(start=str(time0), periods=nc_ini_src0.time.size, freq='1D')
    tref1 = date2num(tref.to_pydatetime(), 'days since 1990-01-01 00:00:00')

    dsaux = nc_roms_grd

    for i in range(nc_ini_src0.time.size):
        output_file = generate_output_file(outfile, tref, i)
        if glob.glob(output_file):
            print(f'{output_file} already saved')
            continue

        # nc_ini_src = nc_ini_src0.isel(time=[i])
        # if glob.glob(outfile % (str(nc_ini_src.time.values[0])[:19])):
        #     print(f'{outfile % (str(nc_ini_src.time.values[0])[:19])} already saved')
        #     continue

        nc_out1 = nc_out0.copy()
        nc_ini_src = nc_ini_src0.isel(time=[i]).rename_dims(rename_coords).rename_vars(rename_vars)
        nc_ini_src.load()
        if horizontal_homog_fields:        

            nc = nc_ini_src.mean(dim=['lon', 'lat'])
            for var in varbs:
                nc_ini_src[var].values[:] = nc[var].values[:, None, None]
            outfile = outfile[:-3] + '_hor_homog.nc'

        zsel = np.delete(np.arange(nc_ini_src.depth.values.size), zdel)
        nc_ini_src_ = nc_ini_src.isel(z=zsel)

        nc_aux1 = interpolate_and_rotate(nc_ini_src_, nc_aux0, dsaux, varbs, dx, dicts, nc_roms_grd)
        copy_boundary_values(nc_aux1, nc_out1, varbs)

        nc_out1 = nc_out1.assign_coords(ocean_time=[tref1[i]])
        nc_out1['ocean_time'].attrs['units'] = 'days since 1990-01-01 00:00:00'
        nc_out1.to_netcdf(output_file)
        nc_out1.close()

if __name__ == '__main__':
    main()