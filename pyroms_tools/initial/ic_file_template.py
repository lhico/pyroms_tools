import pyroms
# from pyroms_tools import remap
from pyroms_tools.utils import utils as ut
from pyroms_tools import utils as ut1

import subprocess
import os
import argparse
import yaml

os.environ["PYROMS_GRIDID_FILE"] = "gridid.txt"

def get_dicts(config_file):
    with open(config_file, 'r') as file:
        dicts = yaml.safe_load(file)
    return dicts['default']

def create_grid_file(dicts):
    ut.create_grid_file(
        dicts['gridid']['gridname'],
        dicts['gridid']['gridname'],
        dicts['grid']['grid'],
        dicts['grid']['N'],
        'roms',
        2,
        dicts['grid']['theta_s'],
        dicts['grid']['theta_b'],
        dicts['grid']['Tcline']
    )

def get_dst_grd(romsgridname):
    return pyroms.grid.get_ROMS_grid(romsgridname)

def create_scalar_fields(dst_grd, starttime):
    for varb in ['zeta', 'temp', 'salt']:
        dfile = f'{dst_grd.name}_{varb}_ic.nc'
        ut1.remap.create_scalar_fields_3D(dst_grd, varb, starttime, dst_file=dfile)

def create_vector_fields(dst_grd, starttime):
    dfileu = f'{dst_grd.name}_u_ic.nc'
    dfilev = f'{dst_grd.name}_v_ic.nc'
    ut1.remap.create_vector_fields(dst_grd, starttime, dst_fileu=dfileu, dst_filev=dfilev)

def merge_files(out_file, ic_file):
    cont = 0
    for istr in ['zeta', 'temp', 'salt', 'u', 'v']:
        if cont == 0:
            command = ('ncks', '--no_abc', '-O', out_file % istr, ic_file)
            cont += 1
        else:
            command = ('ncks', '--no_abc', '-A', out_file % istr, ic_file)
        subprocess.call(command)
        os.remove(out_file % istr)
    print(f'saved to {ic_file}')

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--config', type=str, help='Path to the YAML configuration file')
    args = parser.parse_args()

    dicts = get_dicts(args.config)
    create_grid_file(dicts)
    dst_grd = get_dst_grd(dicts['gridid']['gridname'])
    create_scalar_fields(dst_grd, dicts['ic']['starttime'])
    create_vector_fields(dst_grd, dicts['ic']['starttime'])
    merge_files(f'{dst_grd.name}_%s_ic.nc', dicts['ic']['ic_template'])

if __name__ == "__main__":
    main()