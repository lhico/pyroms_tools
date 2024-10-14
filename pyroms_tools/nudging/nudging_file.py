import argparse
import yaml
from pyroms_tools.utils import utils as ut
import os
import sys

os.environ["PYROMS_GRIDID_FILE"] = "gridid.txt"

def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def create_grid_file(dicts):
    gridfile = dicts['grid']['grid']
    romsgridname = dicts['gridid']['gridname']
    ut.create_grid_file(
        romsgridname,
        romsgridname,
        gridfile,
        dicts['grid']['N'],
        'roms',
        2,
        dicts['grid']['theta_s'],
        dicts['grid']['theta_b'],
        dicts['grid']['Tcline']
    )

def apply_nudging(dicts):
    romsgridname = dicts['gridid']['gridname']
    output = dicts['nudg']['output']
    onudg = ut.nudgcoef(romsgridname)

    east = dicts['nudg']['east']
    west = dicts['nudg']['west']
    north = dicts['nudg']['north']
    south = dicts['nudg']['south']
    tracer_timescales = dicts['nudg']['tracertscale']

    onudg(east, west, north, south, tracer_timescales, foutname=output)

def main():
    parser = argparse.ArgumentParser(description='Process nudging configuration.')
    parser.add_argument('--config', type=str, required=True, help='Path to the config YAML file')
    args = parser.parse_args()

    config = load_config(args.config)
    dicts = config['default']

    create_grid_file(dicts)
    apply_nudging(dicts)

if __name__ == "__main__":
    main()