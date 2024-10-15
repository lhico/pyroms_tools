import os
import sys
import yaml
import logging
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import ndimage
from bathy_smoother import LP_bathy_smoothing, bathy_smoothing
from pyroms import vgrid
from pyroms_toolbox import rx1
from pyroms_tools import utils as ut
import os.path as osp
from argparse import ArgumentParser

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config(config_path, reference):
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config[reference]
    except Exception as e:
        logging.error(f"Error loading configuration: {e}")
        sys.exit(1)

def plot(h2, x, nc0):
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=[10, 3])
    ax = ax.ravel()
    cf = ax[0].contourf(h2.T, levels=np.arange(0, 200, 10), cmap=plt.cm.jet, extend='both')
    ax[1].contourf(nc0.h.T, levels=np.arange(0, 200, 10), cmap=plt.cm.jet, extend='both')
    cf1 = ax[2].contourf((nc0.h - h2).T, levels=np.arange(-200, 210, 10), cmap=plt.cm.RdBu_r)
    ax[3].plot(nc0.h[x])
    ax[3].plot(h2[x])
    ax[0].contour(h2.T, [100, 200, 1000], cmap=plt.cm.jet)
    ax[0].vlines(x, 0, 180)
    ax[1].contour(nc0.h.T, [100, 200, 1000], cmap=plt.cm.jet)
    ax[2].contour(h2.T, [x], colors='r')
    ax[2].contour(nc0.h.T, [x], colors='g')

    for i in range(3):
        ax[i].axis('equal')

    ax[3].set_ylim(4000, 0)

    for i in range(3):
        ax[i].set_xlim(0, 110)
        ax[i].set_ylim(0, 70)

def calculate_rx1(nc, message):
    vgrd = vgrid.s_coordinate_4(nc['h'].values,
                                nc['theta_b'].values,
                                nc['theta_s'].values,
                                nc['Tcline'].values, 
                                nc['s_rho'].size,
                                hraw=nc['hraw'].values)
    z_w = vgrd.z_w.__getitem__(0)
    rmask = nc.mask_rho.values

    logging.info("# -- Haney number -- #\n")
    logging.info("""Steep topography may be a problem in ROMS. Haney number (rx1) is used to assess if topography is too
    steep. 
    rx1~3 is 'safe and conservative'. Values greater than 8-10 are 'insane' 
 (https://www.myroms.org/forum/viewtopic.php?f=14&t=612). Be careful to not oversmooth the bathymetry\n""")
    logging.info(message + '\n')

    rx1values = rx1(z_w, rmask)
    
    return rx1values

def plot_structure(datasetoriginal, datasetnew, fsize=[15, 5], vmin=0, vmax=3000, s=0, vmindiff=-100, vmaxdiff=100):
    datasetoriginal = datasetoriginal.set_coords(['lon_rho', 'lat_rho'])
    datasetnew = datasetnew.set_coords(['lon_rho', 'lat_rho'])
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=fsize, sharex=True, sharey=True)
    datasetoriginal['h'].plot(vmin=vmin, vmax=vmax, x='lon_rho', y='lat_rho', ax=ax[0], cmap=plt.cm.jet)
    datasetnew['h'].plot(vmin=vmin, vmax=vmax, x='lon_rho', y='lat_rho', ax=ax[1], cmap=plt.cm.jet)
    (datasetoriginal['h'] - datasetnew['h']).plot(x='lon_rho', y='lat_rho', ax=ax[2], cmap=plt.cm.RdBu_r, vmin=vmindiff, vmax=vmaxdiff)
    ax[0].set_title('original', loc='right')
    ax[1].set_title('smoothed', loc='right')
    ax[2].set_title('difference', loc='right')
    plt.tight_layout()
    fig.savefig('compare_smoothed_bathymetry_%s.png' % s)



def update_mask(ds):
    ds = ds.copy()
    ds['mask_u'].values = ds['mask_rho'][:,:-1].values * ds['mask_rho'][:,1:].values
    ds['mask_v'].values = ds['mask_rho'][:-1,:].values * ds['mask_rho'][1:,:].values
    ds['mask_psi'].values = ds.mask_u[:-1, :].values * \
                        ds.mask_v[:, :-1].values
    return ds



def main():

    parser = ArgumentParser(description="Process bathymetry data.")
    parser.add_argument('--config', type=str, default='../configs/config.yaml', help='Path to the configuration file')
    parser.add_argument('--reference', type=str, default='default', help='Reference section in the configuration file')
    args = parser.parse_args()

    config_path = args.config
    reference = args.reference

    setup_logging()
    config = load_config(config_path, reference)
    dicts = config['smooth']

    input = dicts['input']
    inputaux = dicts['inputaux']
    output = dicts['output']
    capdepth = dicts['capdepth']
    smooth = dicts['smooth']
    nested_grid = dicts['nested']

    try:
        nc0 = xr.open_dataset(input)
    except Exception as e:
        logging.error(f"Error opening input dataset: {e}")
        sys.exit(1)
    
    nc = nc0.copy()
    nc.load()
    nc = nc.fillna(capdepth)
    nc.mask_rho.values[np.isnan(nc.h.values)] = 0
    nc = update_mask(nc)
    nc = nc.fillna(5)

    rx1in = calculate_rx1(nc, '# -- original rx1 values -- #')

    h = nc.h.values.copy()
    h[h < -capdepth] = -capdepth
    h[h > capdepth] = capdepth

    mask = nc.mask_rho.values.copy()
    h1 = h

    mask0 = mask.copy()
    h1 = bathy_smoothing.smoothing_PlusMinus_rx0(mask, h1, smooth,
            np.gradient(nc.x_rho)[1] * np.gradient(nc.y_rho)[0])

    nc.h.values = h1[0]
    nc.hraw.values[0] = h1[0]
    nc.attrs['smoothing function'] = 'LP_smoothing_rx0'
    nc.attrs['rx0'] = smooth

    rx1out = calculate_rx1(nc, '# -- modified rx1 values -- #')

    nc['rx1'] = nc.h.copy()*0
    nc.rx1.values[:-1, :-1] = rx1out

    nc['smooth_diff'] = nc['hraw'] - nc['h']

    fileout = output
    if osp.exists(fileout):
        os.remove(fileout)
    nc.to_netcdf(fileout)

if __name__ == '__main__':

    main()