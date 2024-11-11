import argparse
import yaml
import sys
import pandas as pd
import xarray as xr
from netCDF4 import date2num
import glob
from pyroms_tools.utils import utils as ut
import numpy as np


def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def process_climatology(config):
    dicts = config['default']
    grid = dicts['grid']['grid']
    source = dicts['clim']['src_file']
    output = dicts['clim']['outfile']
    tstart = dicts['clim']['date'][0]
    tfinal = dicts['clim']['date'][1]
    tfreq = dicts['clim']['date'][2]
    daterange = pd.date_range(start=tstart, end=tfinal, freq=tfreq)

    print('starting loop')
    for datetime in daterange:
        print(datetime)
        dsgrid = xr.open_dataset(grid)
        tref = datetime.to_pydatetime()
        tref1 = date2num(tref, 'days since 1990-01-01 00:00:00')

        output_file = output % str(tref).replace(' ', 'T')

        # if file is already saved continue to the next iteration
        flist = glob.glob(output_file)
        if len(flist) != 0:
            print(f'{flist[0]} already saved')
            continue

        q = ut.interpolationGlorys2Roms(grid, source, tslice=(tstart, tfinal), load=True)
        q.calculate_geostrophy()
        q.set_coords(tref, gtype='rho', time_type=None)
        temp = q.interpolate3d('thetao')
        salt = q.interpolate3d('so')

        q.set_coords(tref, gtype='rho', time_type=None)
        v = q.interpolate3d('vo')
        q.set_coords(tref, gtype='rho', time_type=None)
        u = q.interpolate3d('uo')

        print('interpolation done')

        rot = dsgrid.angle.values
        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v, u)

        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)

        v = q.ds1['vg']
        u = q.ds1['ug']

        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v, u)

        ug = -mag * np.sin(angle + rot)
        vg = mag * np.cos(angle + rot)

        dsgrid = dsgrid.assign_coords(time=[tref])
        dsgrid = dsgrid.assign_coords(temp_time=[tref])
        dsgrid = dsgrid.assign_coords(salt_time=[tref])

        dsgrid['temp'] = (('temp_time', 's_rho', 'eta_rho', 'xi_rho'), [temp])
        dsgrid['salt'] = (('salt_time', 's_rho', 'eta_rho', 'xi_rho'), [salt])
        dsgrid['u'] = (('time', 's_rho', 'eta_u', 'xi_u'), [u1[:, :, :-1]])
        dsgrid['v'] = (('time', 's_rho', 'eta_v', 'xi_v'), [v1[:, :-1, :]])

        dsgrid['ubar'] = (('time', 'eta_u', 'xi_u'), [ug[:, :-1]])
        dsgrid['vbar'] = (('time', 'eta_v', 'xi_v'), [vg[:-1, :]])

        dsgrid['ubar'] = dsgrid['ubar'].bfill('xi_u')
        dsgrid['ubar'] = dsgrid['ubar'].bfill('eta_u')
        dsgrid['vbar'] = dsgrid['vbar'].bfill('xi_v')
        dsgrid['vbar'] = dsgrid['vbar'].bfill('eta_v')

        for v, c in zip(['temp', 'salt', 'v', 'u', 'vbar', 'ubar'], ['rho', 'rho', 'v', 'u', 'v', 'u']):
            barotropic = v in ['ubar', 'vbar']
            dsgrid = ut.nearest_interpolation(dsgrid, v, hgrid=c, barotropic=barotropic)

        dsgrid['time'].attrs = {}
        dsgrid['time'].attrs['long_name'] = 'time'

        dsgrid.to_netcdf(output_file)
        dsgrid.close()
        print(output_file + ' saved')

        del dsgrid, q, ug, vg, u, v

def main():
    parser = argparse.ArgumentParser(description='Process climatology configuration.')
    parser.add_argument('--config', type=str, required=True, help='Path to the config YAML file')
    args = parser.parse_args()

    config = load_config(args.config)
    process_climatology(config)

if __name__ == "__main__":
    main()