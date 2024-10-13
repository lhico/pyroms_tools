from typing import IO
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pyroms
from pyroms_toolbox.CGrid_GLORYS import CGrid_GLORYS
from . import grids as grd

def A2CGrid(grdfile, name='GLORYS_CORAL', area='regional', \
                         xrange=(185,340),
                         yrange=(100, 210),
                         lonoffset=0,
                         latstr='lat',
                         lonstr='lon',
                         depstr='depth',
                         maskstr='mask'):

    nc = xr.open_dataset(grdfile)

    lon = nc[lonstr].values+lonoffset
    lat = nc[latstr].values
    depth = nc[depstr].values
    mask = nc[maskstr].values
    
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
    # depth_bnds[-1] = 6000.
    # # depth_bnds = np.zeros(depth.shape[0]+1)
    # # depth_bnds[:-1] = dep[:,0]
    # # depth_bnds[-1] = dep[-1,1]


    bottom = pyroms.utility.get_bottom(mask_t[::-1], mask_t[0], spval=9999)
    nlev = mask_t.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            if mask_t[0,j,i] == 1:
                h[j,i] =  depth_bnds[int(bottom[j,i])]
    m,l = h.shape

    gridC = CGrid_GLORYS(lon_t, lat_t, lon_u, lat_u, lon_v, lat_v,
        mask_t, mask_u, mask_v, depth_t, depth, h,
        name, xrange, yrange)

    # plt.scatter(gridC.lon_t_vert, gridC.lat_t_vert, marker='.', label='t')
    # plt.scatter(gridC.lon_u_vert, gridC.lat_u_vert, marker='*', label='u')
    # plt.scatter(gridC.lon_v_vert, gridC.lat_v_vert, marker='^', label='v')

    return gridC


class CGrid_GLORYS(object):
    """
    CGrid object for GLORYS
    """

    def __init__(self, lon_t, lat_t, lon_u, lat_u, lon_v, lat_v, mask_t, mask_u, mask_v, depth, depth_bnds, h, name, xrange, yrange):

        self.name = name

        self.xrange = xrange
        self.yrange = yrange

        self.h = h[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t = lon_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_t = lat_t[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_u = lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_u = lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lon_v = lon_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.lat_v = lat_v[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.lon_t_vert = 0.5 * (lon_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_t_vert = 0.5 * (lat_t[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat_t[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.lon_u_vert = 0.5 * (lon_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_u_vert = 0.5 * (lat_u[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat_u[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lon_v_vert = 0.5 * (lon_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lon_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])
        self.lat_v_vert = 0.5 * (lat_v[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1] + \
                               lat_v[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2])

        self.mask_t = mask_t[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_u = mask_u[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        self.mask_v = mask_v[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

        self.z_t = np.tile(depth,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T

        self.z_t_bnds = np.tile(depth_bnds,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T

        ones = np.ones(self.h.shape)
        a1 = lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
             lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a2 = lon_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
             lon_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        a3 = 0.5*(lat_u[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] + \
             lat_u[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1])
        a2 = np.where(a2 > 180*ones, a2 - 360*ones, a2)
        a2 = np.where(a2 < -180*ones, a2 + 360*ones, a2)
        a2 = a2 * np.cos(np.pi/180.*a3)
        self.angle = np.arctan2(a1, a2)
