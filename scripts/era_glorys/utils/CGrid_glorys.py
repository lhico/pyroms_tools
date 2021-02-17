import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pyroms
from pyroms_toolbox.CGrid_GLORYS import CGrid_GLORYS
from . import grids as grd

def A2CGrid(grdfile, name='GLORYS_CORAL', area='regional', \
                         xrange=(185,340), yrange=(100, 210), lonoffset=0):

    nc = xr.open_dataset(grdfile)

    lon = nc.longitude.values+lonoffset
    lat = nc.latitude.values
    depth = nc.depth.values
    mask = nc.mask.values

    gridA = grd.Agrid(lon, lat, depth, mask=mask)
    gridA.A2C()
    # gridA.plot_horizontal_Cgrid(istart=0, iend=30)

    lon_t = gridA.lon_t.copy()
    lat_t = gridA.lat_t.copy()

    print(lon_t.shape)

    lon_u = gridA.lon_u.copy()
    lat_u = gridA.lat_u.copy()
    lon_v = gridA.lon_v.copy()
    lat_v = gridA.lat_v.copy()
    mask_t = gridA.mask_t.copy()
    mask_u = gridA.mask_u.copy()
    mask_v = gridA.mask_v.copy()
    depth_t = gridA.depth_t.copy()
    depth_w = gridA.depth_w.copy()

    depth_bnds = np.zeros(depth_t.shape[0]+1)
    depth_bnds[:-1] = gridA.depth_w[:]
    # depth_bnds[-2] = 6000.
    depth_bnds[-1] = 6000.
    # depth_bnds = np.zeros(depth.shape[0]+1)
    # depth_bnds[:-1] = dep[:,0]
    # depth_bnds[-1] = dep[-1,1]

    # raise ValueError()

    bottom = pyroms.utility.get_bottom(mask_t[::-1], mask_t[0], spval=9999)
    nlev = mask_t.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            if mask_t[0,j,i] == 1:
                h[j,i] = depth_bnds[int(bottom[j,i])]
    m,l = h.shape

    gridC = CGrid_GLORYS(lon_t, lat_t, lon_u, lat_u, lon_v, lat_v,
        mask_t, mask_u, mask_v, depth_t, depth_bnds, h,
        name, xrange, yrange)

    # plt.scatter(gridC.lon_t_vert, gridC.lat_t_vert, marker='.', label='t')
    # plt.scatter(gridC.lon_u_vert, gridC.lat_u_vert, marker='*', label='u')
    # plt.scatter(gridC.lon_v_vert, gridC.lat_v_vert, marker='^', label='v')

    return gridC
