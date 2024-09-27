import matplotlib.pyplot as plt
import os.path as osp
from scipy.interpolate import griddata
import xarray as xr
from scipy import ndimage, signal
import numpy as np
# from dry_package.plot_schemes import maps
import argparse
import cartopy.crs as ccrs
import yaml

def rotate_coords(xm, ym, ang_rot, degrees=False):
    xm_mean = np.mean(xm)
    ym_mean = np.mean(ym)

    xm1 = xm.copy() - xm_mean
    ym1 = ym.copy() - ym_mean

    if degrees:
        ang_rot_rad = ang_rot/180 *np.pi
    ang_rot_rad = np.deg2rad(ang_rot) if degrees else ang_rot

    xrot = xm1 * np.cos(ang_rot_rad) - ym1 * np.sin(ang_rot_rad)
    yrot = xm1 * np.sin(ang_rot_rad) + ym1 * np.cos(ang_rot_rad)

    xrot += xm_mean
    yrot += ym_mean
    return xrot, yrot

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate and rotate grid.')
    parser.add_argument('--config', type=str, default='../configs/config.yaml', help='Path to the config file')
    args = parser.parse_args()

    # -- gets the information from the YAML config file -- #
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)
    dicts = config['default']

    # -- paths -- #
    odir = dicts['output_dir']  # directory where the file will be saved
    
    # -- horizontal grid parameters -- #
    dxdy = dicts['grid']['dxdy']
    x0, x1, y0, y1 = dicts['grid']['WESN']
    xoffset = dicts['grid']['xoffset']
    yoffset = dicts['grid']['yoffset']
    rot = dicts['grid']['rot']

    # -- vertical grid parameters -- #
    N = dicts['grid']['N']
    theta_s = dicts['grid']['theta_s']
    theta_b = dicts['grid']['theta_b']
    Tcline = dicts['grid']['Tcline']

    # -- grid rotation -- #
    x = np.arange(x0, x1, dxdy) + xoffset
    y = np.arange(y0, y1, dxdy) + yoffset

    xm, ym = np.meshgrid(x, y)
    xrot, yrot = rotate_coords(xm, ym, rot)

    # -- define horizontal ROMS grid object -- #
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    # ax = maps.make_overall(extent=[-55, -30, -35, -15])
    ax.scatter(xrot, yrot, marker='.', s=1)
    ax.coastlines(resolution='10m')
    plt.show()