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


if __name__ == '__main__':
    
    # -- gets  the information from the config file -- #
    reference = 'pbs_202109_glorys'
    dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
    dicts = dicts[reference]

    gdir     = dicts['output_dir']
    capdepth = dicts['smooth.capdepth']
    smooth   = dicts['smooth.smooth'] 
    nc0 = xr.open_dataset(osp.join(gdir,'grids/bacia_santos00.nc'))
    nc = nc0.copy()

    # x,y = np.gradient(nc.h)
    # z = (x**2+y**2)**0.5
    h = nc.h.values.copy()
    # smooth = abs(ndimage.gaussian_filter(z/z.max()-1,8))

    h[h<-capdepth] = -capdepth
    h[h>capdepth] = capdepth

    mask = nc.mask_rho.values.copy()
    amp = mask.copy()
    amp[:] = 10000
    # amp[nc0.h<200]=1.
    h1 =  h

    # h1 = bathy_smoothing.smoothing_Positive_rx0(mask, h1, ii)

    # h1 = LP_bathy_smoothing.LP_smoothing_rx0_heuristic(mask,
    #                                         nc0.h.values,ii,
    #                                         -1*np.ones(nc.h.values.shape),
    #                                         amp)

    # h3 = bathy_smoothing.smoothing_Positive_rx0(mask, h2, 0.05)
    h1 = bathy_smoothing.smoothing_PlusMinus_rx0(mask, h1, smooth,
            np.gradient(nc.x_rho)[0]*np.gradient(nc.y_rho)[0])

    # h1 = bathy_smoothing.smoothing_Positive_rx0(mask, h1, ii)


    # h2 = bathy_smoothing.smoothing_PlusMinus_rx0(mask, h1[0], ii,
    #     np.gradient(nc.x_rho)[0]*np.gradient(nc.y_rho)[0])

    nc.h.values = h1[0]
    nc.hraw.values[0] = h1[0]
    nc.attrs['smoothing function'] = 'LP_smoothing_rx0'
    nc.attrs['rx0'] = smooth


    fileout = osp.join(gdir, 'bacia_santos.nc')
    os.system(f'rm {fileout}')
    nc.to_netcdf(fileout)


