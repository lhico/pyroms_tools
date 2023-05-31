from utils import utils as ut
import numpy as np
from bathy_smoother import LP_bathy_smoothing
from bathy_smoother import bathy_smoothing
import xarray as xr
import os
from scipy import ndimage
import matplotlib.pyplot as plt
import os.path as osp
import os
from pyroms import vgrid
from pyroms_toolbox import rx1
import sys


def plot(h2, x):
    fig, ax = plt.subplots(ncols=2,nrows=2,figsize=[10,3])
    ax = ax.ravel()
    cf = ax[0].contourf(h2.T, levels=np.arange(0,200,10), cmap=plt.cm.jet, extend='both')
    ax[1].contourf(nc0.h.T, levels=np.arange(0,200,10), cmap=plt.cm.jet, extend='both')
    cf1 = ax[2].contourf((nc0.h - h2).T, levels=np.arange(-200,210,10), cmap=plt.cm.RdBu_r)
    # ax[2].scatter(s=1, marker='.')
    ax[3].plot(nc0.h[x])
    ax[3].plot(h2[x])
    ax[0].contour(h2.T, [100, 200, 1000], cmap=plt.cm.jet)
    ax[0].vlines(x,0, 180)
    ax[1].contour(nc0.h.T, [100, 200, 1000], cmap=plt.cm.jet)
    ax[2].contour(h2.T, [x], colors='r')
    ax[2].contour(nc0.h.T, [x], colors='g')
    
    for i in range(3):
        ax[i].axis('equal')

    # plt.colorbar(cf1, ax=[ax[2]], orientation='horizontal')
    # plt.colorbar(cf, ax=[ax[0], ax[1]], orientation='horizontal')
    ax[3].set_ylim(4000,0)

    for i in range(3):
        ax[i].set_xlim(0,110)
        ax[i].set_ylim(0,70)


def calculate_rx1(nc,message):

    vgrd = vgrid.s_coordinate_4(nc['h'].values,
                                nc['theta_b'].values,
                                nc['theta_s'].values,
                                nc['Tcline'].values, 
                                nc['s_rho'].size,
                                hraw=nc['hraw'].values)
    z_w = vgrd.z_w.__getitem__(0)
    rmask = nc.mask_rho.values

    print("# -- Haney number -- #\n")
    print("""Steep topography may be a problem in ROMS. Haney number (rx1) is used to asses if topography is too
    steep. 
    rx1~3 is 'safe and conservative'. Values greater than 8-10 are 'insane' 
 (https://www.myroms.org/forum/viewtopic.php?f=14&t=612). Be careful to not oversmooth the bathymetry\n""")
    print(message + '\n')

    rx1values = rx1(z_w,rmask)
    
    return rx1values


def plot_structure(datasetoriginal, datasetnew, fsize=[15,5], vmin=0, vmax=3000, s=0,vmindiff=-100,vmaxdiff=100):
    datasetoriginal = datasetoriginal.set_coords(['lon_rho', 'lat_rho'])
    datasetnew = datasetnew.set_coords(['lon_rho', 'lat_rho'])
    fig, ax = plt.subplots(nrows=1,ncols=3,figsize=fsize,sharex=True,sharey=True)
    datasetoriginal['h'].plot(vmin=vmin, vmax=vmax,x='lon_rho', y='lat_rho', ax=ax[0], cmap=plt.cm.jet)
    datasetnew['h'].plot(vmin=vmin, vmax=vmax,x='lon_rho', y='lat_rho', ax=ax[1], cmap=plt.cm.jet)
    (datasetoriginal['h']-datasetnew['h']).plot(x='lon_rho', y='lat_rho', ax=ax[2],cmap=plt.cm.RdBu_r, vmin=vmindiff, vmax=vmaxdiff)
    ax[0].set_title('original', loc='right')
    ax[1].set_title('smoothed', loc='right')
    ax[2].set_title('difference', loc='right')
    plt.tight_layout()
    fig.savefig('compare_smoothed_bathymetry_%s.png' % s)

if __name__ == '__main__':
    
    nested_grid = True
    # -- gets  the information from the config file -- #
    reference = sys.argv[1] if len(sys.argv)==2  else 'swatl_2022_deep4_nested'

    dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
    dicts = dicts[reference]

    input       = dicts['smooth.input']
    inputaux       = dicts['smooth.inputaux']
    output        = dicts['smooth.output']
    capdepth    = dicts['smooth.capdepth']
    smooth      = dicts['smooth.smooth']
    nested_grid = dicts['smooth.nested']

    nc0 = xr.open_dataset(input)
    
    # matlab function does not copy some values from original file
    # to nested file (Matlab R2021b Update 4 (9.11.0.2022996) 64-bit (glnxa64))
    if nested_grid:
        ncaux = xr.open_dataset(inputaux) 
        nc0.load()
        for i in ['theta_s', 'theta_b','hc', 'Cs_r', 'Cs_w', 'Tcline']:
            nc0[i].values = ncaux[i].values
        nc0.assign_coords(s_rho=ncaux.s_rho)
    

    nc = nc0.copy()
    nc = nc.fillna(capdepth)
    nc.mask_rho.values[np.isnan(nc.h.values)] = 0
    nc = ut.update_mask(nc)
    nc = nc.fillna(5)


    rx1in = calculate_rx1(nc, '# -- original rx1 values -- #')

    h = nc.h.values.copy()
    # smooth = abs(ndimage.gaussian_filter(z/z.max()-1,8))

    h[h<-capdepth] = -capdepth
    h[h>capdepth] = capdepth

    mask = nc.mask_rho.values.copy()

    h1 =  h

    # h1 = bathy_smoothing.smoothing_Positive_rx0(mask, h1, ii)

    mask0 = mask.copy()

    h1 = bathy_smoothing.smoothing_PlusMinus_rx0(mask, h1, smooth,
            np.gradient(nc.x_rho)[0]*np.gradient(nc.y_rho)[0])


    nc.h.values = h1[0]
    nc.hraw.values[0] = h1[0]
    nc.attrs['smoothing function'] = 'LP_smoothing_rx0'
    nc.attrs['rx0'] = smooth

    rx1out = calculate_rx1(nc, '# -- modified rx1 values -- #')

    nc['rx1'] = nc.h.copy()
    nc.rx1.values[:-1,:-1] = rx1out 

    # fileout = osp.join(gdir, 'bacia_santos.nc')
    fileout = outdir
    os.system(f'rm {fileout}')
    nc.to_netcdf(fileout)
