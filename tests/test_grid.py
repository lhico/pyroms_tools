import pathlib
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import subprocess
import pytest
import pyroms
from pyroms_tools.grid.make_grid_ import (load_config_from_yaml, rotate_coords, 
                                          interpolate_bathymetry, hgrid, h_bathymetry)

@pytest.fixture
def expected_grid():
    basedir = pathlib.Path(__file__).parent.resolve()
    return xr.open_dataset(f'{basedir}/expected/expected_grid.nc')

@pytest.fixture
def test_grid():
    """
        function to create a numerical domain with real bathymetry based on a small piece of 
        the GEBCO2019 dataset
    """
    basedir = pathlib.Path(__file__).parent.resolve()

    dicts = load_config_from_yaml(f"{basedir}/grid_config_test.yml",
                                  "default")

    # Extract parameters from the configuration
    bfile = dicts['bathy_file']
    dxdy = dicts['grid']['dxdy']
    x0, x1, y0, y1 = dicts['grid']['WESN']
    xoffset = dicts['grid']['xoffset']
    yoffset = dicts['grid']['yoffset']
    rot = dicts['grid']['rot']
    n = dicts['grid']['N']
    theta_s = dicts['grid']['theta_s']
    theta_b = dicts['grid']['theta_b']
    tcline = dicts['grid']['Tcline']
    gridout = dicts['grid']['grid']
    grd_name = dicts['grid']['grid_name']

    # Generate grid coordinates
    x = np.arange(x0, x1, dxdy) + xoffset
    y = np.arange(y0, y1, dxdy) + yoffset
    xm, ym = np.meshgrid(x, y)
    xrot, yrot = rotate_coords(xm, ym, rot)

    # Create horizontal grid and interpolate bathymetry
    hgrd, _ = hgrid(xrot, yrot)
    topo = interpolate_bathymetry(bfile, hgrd)
    lon = xrot
    lat = yrot
    h, hraw = h_bathymetry(topo, lon, lat, hgrd)

    # Close all plots
    plt.close('all')

    # Create vertical grid and ROMS grid
    vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, tcline, n, hraw=hraw)
    grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

    # Write the grid to a file
    pyroms.grid.write_ROMS_grid(grd, filename=gridout)

    return xr.open_dataset(f'{basedir}/test_grid.nc')

def testing_fixture_grids(expected_grid, test_grid):   
    # assert both
    assert expected_grid.equals(test_grid), "Numerical grids are not equal."

    basedir = pathlib.Path(__file__).parent.resolve()
    subprocess.call(f'rm {basedir}/test_grid.nc', shell=True)
