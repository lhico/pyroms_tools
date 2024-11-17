import pathlib
import subprocess
import xarray as xr

import pytest

from pyroms_tools.initial.ic_file_template import (get_dst_grd, create_scalar_fields, 
                                                   create_vector_fields)
from pyroms_tools.utils.utils import create_grid_file

from test_utils import load_yaml_with_includes

@pytest.fixture
def expected_initial_template():
    basedir = pathlib.Path(__file__).parent.resolve()
    return xr.open_dataset(f'{basedir}/expected/expected_initial_template.nc')

@pytest.fixture
def test_initial_template():
    basedir = pathlib.Path(__file__).parent.resolve()

    dicts = load_yaml_with_includes(f"{basedir}/initial_config_test.yml")

    create_grid_file(
        dicts['grid']['grid']['gridname'],
        dicts['grid']['grid']['gridname'],
        dicts['grid']['grid']['grid'].replace('/test_grid.nc', '/expected/expected_grid.nc'),
        dicts['grid']['grid']['N'],
        'roms',
        2,
        dicts['grid']['grid']['theta_s'],
        dicts['grid']['grid']['theta_b'],
        dicts['grid']['grid']['Tcline']
    )

    dst_grd = get_dst_grd(dicts['grid']['grid']['gridname'])
    create_scalar_fields(dst_grd, dicts['ic']['starttime'])
    create_vector_fields(dst_grd, dicts['ic']['starttime'])

    ds = xr.open_mfdataset(f'{dst_grd.name}_*_ic.nc')

    subprocess.call(f'rm {basedir}/../{dst_grd.name}_*_ic.nc', shell=True)
    subprocess.call(f'rm {basedir}/../gridid.txt', shell=True)

    return ds

def testing_fixture_initial_templates(expected_initial_template, test_initial_template):   
    # assert both
    assert expected_initial_template.equals(test_initial_template), "Initial condition template are not equal."
