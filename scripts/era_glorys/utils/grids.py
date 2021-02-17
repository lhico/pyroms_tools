import numpy as np
import pyroms
import matplotlib.pyplot as plt
from pyroms_toolbox.CGrid_GLORYS import CGrid_GLORYS

# my bugged CGRID
# class CGrid(object):
#     """
#     CGrid object
#     """
#
#     def __init__(self, lon_t, lat_t, lon_u, lat_u, lon_v, lat_v,
#                  mask_t, mask_u, mask_v, depth, depth_bnds, h,
#                  name, xrange, yrange):
#
#
#         self.name = name
#
#         self.xrange = xrange
#         self.yrange = yrange
#
#         # -- cutting original grid -- #
#         self.h = self._get_ranges(h)
#
#         self.lon_t_aux = self._get_ranges(lon_t)
#         self.lat_t_aux = self._get_ranges(lat_t)
#
#         self.lon_u_aux = self._get_ranges(lon_u)
#         self.lat_u_aux = self._get_ranges(lat_u)
#
#         self.lon_v_aux = self._get_ranges(lon_v)
#         self.lat_v_aux = self._get_ranges(lat_v)
#
#         # -- getting vert for t, u, v --#
#         self.lon_t_vert = self._get_vert1(lon_t)
#         self.lat_t_vert = self._get_vert1(lat_t)
#
#         self.lon_u_vert = self._get_vert1(lon_u)
#         self.lat_u_vert = self._get_vert1(lat_u)
#
#         self.lon_v_vert = self._get_vert1(lon_v)
#         self.lat_v_vert = self._get_vert1(lat_v)
#
#         #-- getting mask ranges --#
#         self.mask_t = self._get_ranges(mask_t)
#         self.mask_u = self._get_ranges(mask_u)
#         self.mask_v = self._get_ranges(mask_v)
#
#         self.z_t = np.tile(depth,
#             (self.mask_t.shape[2],self.mask_t.shape[1],1)).T
#
#         self.z_t_bnds = np.tile(
#             depth_bnds,(self.mask_t.shape[2],self.mask_t.shape[1],1)).T
#
#         #-- setting up angles in the grid --#
#         ones = np.ones(self.h.shape)
#         a1 = self._get_vert2(lat_u)
#         a2 = self._get_vert2(lon_u)
#         a3 = 0.5*self._get_vert2(lat_u, subtract=True)
#
#         a2 = np.where(a2 > 180*ones, a2 - 360*ones, a2)
#         a2 = np.where(a2 < -180*ones, a2 + 360*ones, a2)
#         a2 = a2 * np.cos(np.pi/180.*a3)
#         self.angle = np.arctan2(a1, a2)
#
#
#     def _get_vert1(self, varb):
#         xrange = self.xrange
#         yrange = self.yrange
#         varb1 = varb.copy()
#         varbTVert = varb1[yrange[0]-1:yrange[1]+1, xrange[0]-1:xrange[1]+1]
#         varbTVert += varb1[yrange[0]:yrange[1]+2, xrange[0]:xrange[1]+2]
#         varbTVert *= 0.5
#         return varbTVert
#
#     def _get_vert2(self, varb, subtract=False):
#         xrange = self.xrange
#         yrange = self.yrange
#         varb1 = varb.copy()
#
#         if subtract:
#             varbTVert = varb1[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] - \
#                         varb1[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
#         else:
#             varbTVert = varb1[yrange[0]:yrange[1]+1, xrange[0]+1:xrange[1]+2] + \
#                        varb1[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
#         return varbTVert
#
#     def _get_ranges(self, varb):
#         xrange = self.xrange
#         yrange = self.yrange
#         varb1 = varb.copy()
#
#         if varb.ndim == 2:
#             return varb[yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
#         elif varb.ndim == 3:
#             return varb[:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
#         else:
#             raise IOError('input can have 2 or 3 dimensions')
#
#
#     def plot_horizontal_Cgrid(self, fig=None, ax=None, figkw={}, axkw={}):
#
#         if fig is None:
#             fig = plt.figure(**figkw)
#         if ax is None:
#             ax = fig.add_subplot(**axkw)
#
#         lont = self.lon_t_vert
#         latt = self.lat_t_vert
#
#         lonu = self.lon_u_vert
#         latu = self.lat_u_vert
#
#         lonv = self.lon_v_vert
#         latv = self.lat_v_vert
#
#         ax.scatter(lont, latt, marker='.', label='t')
#         ax.scatter(lonu, latu, marker='*', label='u')
#         ax.scatter(lonv, latv, marker='^', label='v')
#         plt.legend()
#
#         return ax


class Agrid(object):
    def __init__(self, lon, lat, depth, mask=None):
        if lat.ndim + lon.ndim + depth.ndim > 3:
            raise IOError('lat, lon and depth must be unidimensional.' \
            'Curvilinear grids are not implemented yet')

        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.mask = mask

    def A2C(self):
        lon = self.lon
        lat = self.lat

        depth = self.depth
        # depth = self.depth[:-1]  # quickfix

        lont, latt = np.meshgrid(lon, lat)

        self.depth_t = depth

        self.lon_t = lont
        self.lat_t = latt

        self.lon_u = lont - np.gradient(lont, axis=1)/2
        self.lat_u = latt

        self.lon_v = lont
        self.lat_v = latt - np.gradient(latt, axis=0)/2

        self.depth_w = depth + np.gradient(depth, axis=0)/2

        if self.mask is not None:
            mask = self.mask
            self.mask_t = mask
            self.mask_u = mask
            self.mask_v = mask

    def plot_horizontal_Cgrid(self, istart=0, iend=-1, dx=1, dy=1,
        fig=None, ax=None, figkw={}, axkw={}):
        self.A2C()

        if fig is None:
            fig = plt.figure(**figkw)
        if ax is None:
            ax = fig.add_subplot(**axkw)

        lont = self.lon_t
        latt = self.lat_t

        lonu = self.lon_u
        latu = self.lat_u

        lonv = self.lon_v
        latv = self.lat_v

        ax.scatter(lont[istart:iend:dy, istart:iend:dx], latt[istart:iend:dy, istart:iend:dx], marker='.', label='t')
        ax.scatter(lonu[istart:iend:dy, istart:iend:dx], latu[istart:iend:dy, istart:iend:dx], marker='*', label='u')
        ax.scatter(lonv[istart:iend:dy, istart:iend:dx], latv[istart:iend:dy, istart:iend:dx], marker='^', label='v')
        plt.legend()

        return ax



if __name__ == '__main__':

    name = 'test'
    xrange=(2,18)
    yrange=(2,18)

    x = np.arange(0,0.1, 0.001).astype('float32')
    y = np.arange(0,0.1, 0.001).astype('float32')
    z = np.arange(20).astype('float32')
    mask = np.ones([20,x.size,y.size])

    gridA = Agrid(x, y, z, mask=mask)
    gridA.A2C()
    # gridA.plot_horizontal_Cgrid()

    lon_t = gridA.lon_t.copy()
    lat_t = gridA.lat_t.copy()
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
    depth_bnds[-1] = 2.

    bottom = pyroms.utility.get_bottom(mask_t[::-1], mask_t[0], spval=9999)
    nlev = mask_t.shape[0]
    bottom = (nlev-1) - bottom
    h = np.zeros(mask_t[0,:].shape)
    for i in range(mask_t[0,:].shape[1]):
        for j in range(mask_t[0,:].shape[0]):
            # if mask_t[0,j,i] == 1:
            h[j,i] = depth_bnds[int(bottom[j,i])]

    gridC_GLORYS = CGrid_GLORYS(lon_t, lat_t, lon_u, lat_u, lon_v, lat_v,
        mask_t, mask_u, mask_v, depth_t, depth_bnds, h,
        name, xrange, yrange)

    plt.scatter(gridC_GLORYS.lon_t_vert, gridC_GLORYS.lon_t_vert, marker='.', label='t')
    plt.scatter(gridC_GLORYS.lon_u_vert, gridC_GLORYS.lon_u_vert, marker='*', label='u')
    plt.scatter(gridC_GLORYS.lon_v_vert, gridC_GLORYS.lon_v_vert, marker='^', label='v')
