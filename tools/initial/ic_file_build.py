from typing import Tuple, Dict, Any
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
from scipy import interpolate
from utils import utils as ut
import os, sys
from scipy.spatial import cKDTree
from netCDF4 import date2num, num2date
import pandas as pd
import datetime as dtt
import glob
import logging
import argparse
import yaml

# Configure logging
logging.basicConfig(level=logging.INFO)

def compute_depth_layers(ds: xr.Dataset, hmin: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """Compute depths of ROMS vertical levels (Vtransform = 2)."""
    S_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
    S_w = (ds.hc * ds.s_w + ds.Cs_w * ds.h) / (ds.hc + ds.h)
    z_rho = ds.h * S_rho
    z_w = ds.h * S_w
    return z_rho, z_w

def setup_nc_roms_grd(nc_roms_grd: xr.Dataset, nc0: xr.Dataset) -> xr.Dataset:
    """Setup ROMS grid with coordinates and values from another dataset."""
    nc_roms_grd = nc_roms_grd.assign_coords(s_rho=nc0.s_rho)
    nc_roms_grd = nc_roms_grd.assign_coords(s_w=nc0.s_w)
    nc_roms_grd.Cs_r.values = nc0.Cs_r.values
    nc_roms_grd.hc.values = nc0.hc.values
    nc_roms_grd.Cs_w.values = nc0.Cs_w.values
    return nc_roms_grd

def interpolate_horizontal(source_grid: xr.Dataset, target_grid: xr.Dataset) -> xr.Dataset:
    """Interpolate data from source grid to target grid using bilinear method."""
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear', extrap_method='nearest_s2d')
    return regridder(source_grid)

def interpolation(fpath: str, nc_roms_grd: xr.Dataset, source_grid: xr.Dataset, target_grid: xr.Dataset, gridtype: str = 'rho') -> np.ndarray:
    """Interpolate data vertically and horizontally."""
    target_grid = target_grid.copy()
    coords_rename = get_coords_rename()

    nc0 = xr.open_dataset(fpath)
    nc_roms_grd = setup_nc_roms_grd(nc_roms_grd, nc0)
    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])
    z, _ = compute_depth_layers(nc_roms_grd)
    z = z.transpose(*('s_rho', 'eta_rho', 'xi_rho'), transpose_coords=False).values

    if gridtype == 'u':
        z = z[:, :, :-1]
    elif gridtype == 'v':
        z = z[:, :-1, :]

    interpolated = interpolate_horizontal(source_grid, target_grid)
    interpvarb = np.zeros(z.shape)
    mask = nc_roms_grd[f'mask_{gridtype}'].values
    ind = np.where(mask != 0)

    for j, i in zip(ind[0], ind[1]):
        logging.info(f'Interpolating: {j}, {i}')
        f = interpolate.interp1d(-interpolated.depth.values,
                                 interpolated[:, j, i].values,
                                 bounds_error=False,
                                 fill_value='extrapolate',
                                 kind='slinear')
        interpvarb[:, j, i] = f(z[:, j, i])
    return interpvarb

def interpolation2d(fpath: str, nc_roms_grd: xr.Dataset, source_grid: xr.Dataset, target_grid: xr.Dataset, gridtype: str = 'rho') -> np.ndarray:
    """Interpolate data horizontally."""
    target_grid = target_grid.copy()
    coords_rename = get_coords_rename()

    nc0 = xr.open_dataset(fpath)
    nc_roms_grd = setup_nc_roms_grd(nc_roms_grd, nc0)
    target_grid = target_grid.rename(coords_rename['roms2data'][gridtype])
    z, _ = compute_depth_layers(nc_roms_grd)
    z = z.transpose(*('s_rho', 'eta_rho', 'xi_rho'), transpose_coords=False).values
    interpolated = interpolate_horizontal(source_grid, target_grid)
    return interpolated.values

def extrapolation_nearest(x: np.ndarray, y: np.ndarray, var: np.ndarray, maskvalue: Any = None) -> np.ndarray:
    """Nearest-neighbor extrapolation with cKDTree method."""
    varx = var.copy()
    if maskvalue is not None:
        varx[var == maskvalue] = np.nan
    varz = varx.copy()
    if x.ndim == 1:
        x, y = np.meshgrid(x, y)
    n = var.shape[0]
    for k in range(var.shape[0]):
        logging.info(f'Nearest extrapolation: {k/n*100:0.2f}%')
        idxnan = np.where(np.isnan(varx[k]))
        idx = np.where(~np.isnan(varx[k]))
        wet = np.zeros((len(idx[0]), 2)).astype(int)
        dry = np.zeros((len(idxnan[0]), 2)).astype(int)
        wet[:, 0] = idx[0].astype(int)
        wet[:, 1] = idx[1].astype(int)
        dry[:, 0] = idxnan[0].astype(int)
        dry[:, 1] = idxnan[1].astype(int)
        xwet = x[wet[:, 0], wet[:, 1]]
        ywet = y[wet[:, 0], wet[:, 1]]
        xdry = x[dry[:, 0], dry[:, 1]]
        ydry = y[dry[:, 0], dry[:, 1]]
        xy = np.array([xdry, ydry])
        tree = cKDTree(np.c_[xwet, ywet])
        _, ii = tree.query(xy.T, k=1)
        if wet.shape[0] == 0:
            pass
        else:
            varz[k, dry[:, 0], dry[:, 1]] = varx[k, wet[ii, 0], wet[ii, 1]]
    return varz
def regrid(source_grid, target_grid):
    regridder = xe.Regridder(source_grid, target_grid, 'bilinear')
    return regridder(source_grid)

def rotate_vector_field(nc: xr.Dataset) -> xr.Dataset:
    """Rotate vector fields in ROMS IC file."""
    nc_out = nc.copy()
    nc_out.load()
    
    target_grid = nc_out.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
    
    interp_u = regrid(nc_out.rename({'lon_u': 'lon', 'lat_u': 'lat'}), target_grid)
    interp_v = regrid(nc_out.rename({'lon_v': 'lon', 'lat_v': 'lat'}), target_grid)
    
    rot = nc_out.angle
    for var in [['u', 'v'], ['ubar', 'vbar']]:
        u = interp_u[var[0]]
        v = interp_v[var[1]]
        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v, u)
        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)
        interp_u[var[0]] = u1
        interp_v[var[1]] = v1
    
    regrid_u = regrid(interp_u, nc_out.rename({'lon_u': 'lon', 'lat_u': 'lat'}))
    regrid_v = regrid(interp_v, nc_out.rename({'lon_v': 'lon', 'lat_v': 'lat'}))
    
    for var in [['u', 'v'], ['ubar', 'vbar']]:
        nc_out[var[0]] = regrid_u[var[0]]
        nc_out[var[1]] = regrid_v[var[1]]
    
    return nc_out


def get_coords_rename() -> Dict[str, Dict[str, Dict[str, str]]]:
    """Get coordinate renaming mappings."""
    return {
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



def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


def main():
    parser = argparse.ArgumentParser(description="Process initial conditions for ROMS.")
    parser.add_argument('--config', type=str, help='Path to the configuration YAML file.')
    parser.add_argument('--reference', type=str, default='ceresIV_2.012', help='Reference domain from the config file.')
    args = parser.parse_args()

    config = read_config(args.config)
    dicts = config[args.reference]

    outfile = dicts['ic']['ic_template']
    ic_start = dicts['ic']['starttime']
    rename_coords = dicts.get('rename_dims', {})
    rename_vars = dicts.get('rename_vars', {})
    varbs = dicts['ic']['interp_varbs']
    map_varbs = dicts['ic']['map_varbs']
    invert_depth = dicts.get('invert', [])
    zdel = dicts.get('delete_idepths', [])
    horizontal_homog_fields = dicts.get('ic.hor_homog', False)

    nc_roms_grd = xr.open_dataset(dicts['grid']['grid'])
    nc_ini_src = xr.open_mfdataset(dicts['ic']['source_file'], decode_times=False, chunks={'time': 1})
    nc_out0 = xr.open_dataset(dicts['ic']['ic_file'])
    nc_out = nc_out0.copy()

    ic_start = dtt.datetime.strptime(ic_start, '%Y-%m-%dT%H:%M:%S')
    tref = pd.date_range(start=str(ic_start), periods=1, freq='1H')
    tref1 = date2num(tref.to_pydatetime(), 'days since 1990-01-01 00:00:00')

    if len(nc_ini_src.dims) == 4:
        if ic_start is None:
            nc_ini_src = nc_ini_src.isel(time=[0])
        else:
            nc_ini_src = nc_ini_src.sel(time=[tref1], method='nearest')

    outfile = ic_start.strftime(outfile)

    if glob.glob(outfile):
        print(f'{outfile} already saved')
        exit()

    ds_out = nc_roms_grd

    nc_ini_src = nc_ini_src.rename_dims(rename_coords).rename_vars(rename_vars)
    print(nc_ini_src.time.values)

    if horizontal_homog_fields:
        nc = nc_ini_src.mean(dim=['lon', 'lat'])
        for var in varbs:
            print(var)
            if nc_ini_src[var].ndim == 3:
                if var in ['salt', 'temp']:
                    raise IOError(f'{var} should be 4D')
                nc_ini_src[var].values[:] = 0
            elif nc_ini_src[var].ndim == 4:
                nc_ini_src[var].values[0, :] = nc[var].values[0, :, None, None]
            nc_ini_src['u'].values[:] = 0
            nc_ini_src['v'].values[:] = 0
            nc_ini_src['zeta'].values[:] = 0

        outfile = outfile[:-3] + '_hor_homog.nc'

    zsel = np.arange(nc_ini_src.depth.values.size)
    zsel = np.delete(zsel, zdel)
    nc_ini_src = nc_ini_src.isel(depth=zsel)

    for varb in varbs:
        print(f'\n# -- Interpolating {varb} --#\n')

        gtype = 'rho'
        if varb == 'uo':
            gtype = 'u'
        elif varb == 'vo':
            gtype = 'v'

        print(f'Grid type: {gtype}')
        print(f'{varb} shape: {nc_ini_src[varb].values.ndim}')

        if "time" not in nc_ini_src[varb].coords:
            raise IOError("Time should be a coordinate. If there is a time coordinate with a different name, please rename it to 'time'")

        nc_ini_src[varb].load()
        if nc_ini_src[varb].values.ndim == 4:
            nc_ini_src[varb].values[0] = extrapolation_nearest(nc_ini_src.longitude.values,
                                                               nc_ini_src.latitude.values,
                                                               nc_ini_src[varb].values[0,])
            var = interpolation(dicts['grid']['grid'], nc_roms_grd, nc_ini_src[varb][0], ds_out, gridtype=gtype)
        elif nc_ini_src[varb].values.ndim == 3:
            nc_ini_src[varb].values = extrapolation_nearest(nc_ini_src.longitude.values,
                                                            nc_ini_src.latitude.values,
                                                            nc_ini_src[varb].values)
            var = interpolation2d(dicts['grid']['grid'], nc_roms_grd, nc_ini_src[varb], ds_out, gridtype=gtype).squeeze()
        else:
            raise IndexError('Array should be either 3D or 4D')

        nc_out[map_varbs[varb]].values = [var]

        if nc_out[map_varbs[varb]].ndim == 4:
            nc_out[map_varbs[varb]] = nc_out[map_varbs[varb]].bfill(dim='s_rho').load()

    nc_out1 = rotate_vector_field(nc_out)

    nc_out1.ubar.values[:] = ((nc_out1.u * nc_out1.s_rho).sum(dim='s_rho') / nc_out1.s_rho.sum()).values
    nc_out1.vbar.values[:] = ((nc_out1.v * nc_out1.s_rho).sum(dim='s_rho') / nc_out1.s_rho.sum()).values
    nc_out1 = nc_out1.assign_coords(ocean_time=tref1)
    nc_out1['ocean_time'].attrs['units'] = 'days since 1990-01-01 00:00:00'

    os.system(f'rm {outfile}')
    nc_out1.to_netcdf(outfile)

if __name__ == "__main__":
    main()