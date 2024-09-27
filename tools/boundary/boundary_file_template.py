import pyroms
from utils import CGrid_glorys
from utils import remap_bdry
import subprocess
import os
import os.path as osp
import pandas as pd
import sys
import yaml
from utils import utils as ut
import argparse

os.environ["PYROMS_GRIDID_FILE"] = "gridid.txt"

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some parameters.')
    parser.add_argument('--config', type=str, help='Path to the YAML configuration file')
    return parser.parse_args()

def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

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

def get_timerange(startTime):
    return pd.date_range(startTime, startTime, freq='1D')

def interpolate_vector_fields(dst_grd, itstr):
    dfileu = f'{dst_grd.name}_u_bdry.nc'
    dfilev = f'{dst_grd.name}_v_bdry.nc'
    remap_bdry.make_bdry_uv_file(dst_grd, itstr, dst_fileu=dfileu, dst_filev=dfilev)

def interpolate_scalar_fields(dst_grd, itstr):
    for varb in ['so', 'thetao', 'zos']:
        dfile = f'{dst_grd.name}_{varb}_bdry.nc'
        remap_bdry.make_bdry_rho_file(dst_grd, varb, itstr, dst_file=dfile)

def merge_files_and_cleanup(output_file, dst_grd):
    out_file = f'{dst_grd.name}_%s_bdry.nc'
    cont = 0
    for istr in ['zos', 'so', 'thetao', 'u', 'v']:
        if cont == 0:
            command = ('ncks', '--no_abc', '-O', out_file % istr, output_file)
            cont += 1
        else:
            command = ('ncks', '--no_abc', '-A', out_file % istr, output_file)
        subprocess.call(command)
        os.remove(out_file % istr)

def main():
    args = parse_arguments()
    config = load_config(args.config)
    
    romsgridname = config['ceresIV_2.012']['gridid']['gridname']
    startTime = config['ceresIV_2.012']['bndry']['startTime']  # This should be updated based on your actual config
    output_file = config['ceresIV_2.012']['bndry']['bdry_template']
    gridfile = config['ceresIV_2.012']['grid']['grid']
    
    create_grid_file(config['ceresIV_2.012'])
    
    timerange = get_timerange(startTime)
    dst_grd = pyroms.grid.get_ROMS_grid(romsgridname)
    
    it = timerange[0]
    itstr = str(it).replace(' ', 'T')
    
    interpolate_vector_fields(dst_grd, itstr)
    interpolate_scalar_fields(dst_grd, itstr)
    
    merge_files_and_cleanup(output_file, dst_grd)

if __name__ == "__main__":
    main()