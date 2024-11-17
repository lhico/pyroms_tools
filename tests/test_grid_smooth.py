import pathlib, os
import os.path as osp
import subprocess

import numpy as np
import xarray as xr

import pytest

from bathy_smoother import bathy_smoothing

from pyroms_tools.grid.make_grid_smooth import (
    load_config,
    calculate_rx1,
    update_mask
)

@pytest.fixture
def expected_smooth():
    basedir = pathlib.Path(__file__).parent.resolve()
    return xr.open_dataset(f'{basedir}/expected/expected_grid_smooth.nc')

@pytest.fixture
def test_smooth():
    basedir = pathlib.Path(__file__).parent.resolve()

    dicts = load_config(f"{basedir}/grid_config_test.yml",
                                  "default")

    inputgrid = dicts['smooth']['input']
    output = dicts['smooth']['output']
    capdepth = dicts['smooth']['capdepth']
    smooth = dicts['smooth']['smooth']
    
    nc0 = xr.open_dataset(inputgrid)

    nc = nc0.copy()
    nc.load()
    nc = nc.fillna(capdepth)
    nc.mask_rho.values[np.isnan(nc.h.values)] = 0
    nc = update_mask(nc)
    nc = nc.fillna(5)

    h = nc.h.values.copy()
    h0 = nc.h.values.copy()
    
    h[h < -capdepth] = -capdepth
    h[h > capdepth] = capdepth

    mask = nc.mask_rho.values.copy()
    h1 = h

    h1 = bathy_smoothing.smoothing_PlusMinus_rx0(mask, h1, smooth,
            np.gradient(nc.x_rho)[1] * np.gradient(nc.y_rho)[0])

    nc.h.values = h1[0]
    nc.hraw.values[0] = h1[0]
    nc.attrs['smoothing function'] = 'LP_smoothing_rx0'
    nc.attrs['rx0'] = smooth

    rx1out = calculate_rx1(nc, '# -- modified rx1 values -- #')

    nc['rx1'] = nc.h.copy()*0
    nc.rx1.values[:-1, :-1] = rx1out

    nc['smooth_diff'] = h0 - nc['h']

    fileout = output
    if osp.exists(fileout):
        os.remove(fileout)
    nc.to_netcdf(fileout)

    return xr.open_dataset(fileout)

def testing_fixture_smooth(expected_smooth, test_smooth):   
    # assert both
    assert expected_smooth.equals(test_smooth), "Numerical smoothed grids are not equal."

    basedir = pathlib.Path(__file__).parent.resolve()
    subprocess.call(f'rm {basedir}/test_grid_smooth.nc', shell=True)
