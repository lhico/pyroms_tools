import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy as sp

auxipath  = '../../data/aux_data'
romspath = '../../data/roms_files'

varbdict = {
    'roms': {
        'depth': 'h'
    },
    'merc': {
        'depth': 'deptho'
    }
}

ncmerc = xr.open_dataset(auxipath + '/glo_mask_bathy2.nc')
ncroms = xr.open_dataset(romspath + '/sbb_grid_roms_100m.nc')

lonmerc = ncmerc['longitude']+ 360
latmerc = ncmerc['latitude']
lonroms = ncroms['lon_rho']
latroms = ncroms['lat_rho']
depmerc = ncmerc[varbdict['merc']['depth']]
deproms = ncroms[varbdict['roms']['depth']]


cvalues = [0, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 2000, 4000]
cvalues = np.arange(2000,5000, 50)

sigma_y = 0.75
sigma_x = 0.75

# interpolate
lon = lonroms.values
lat = latroms.values
lonref = lonmerc.values
latref = latmerc.values
xm, ym = np.meshgrid(lonref, latref)


# interpolation.tree.query(zip(xm.ravel(), ym.ravel()))
bla = interpolate.griddata(list(zip(xm.ravel(), ym.ravel())),
                     depmerc.values.ravel(),
                     list(zip(lon.ravel(), lat.ravel())))

bla = bla.reshape(lon.shape)
bla[np.isnan(bla)] =5
# plt.contourf(bla - deproms.values, np.arange(0,100))

yi_weights = np.cosh(np.arange(-76,77))**0.25
xi_weights = np.cosh(np.arange(-111,112))**0.25

yi_weights = np.delete(yi_weights, int(yi_weights.size/2))
xi_weights = np.delete(xi_weights, int(xi_weights.size/2))

yi_weights = yi_weights/yi_weights.max()
xi_weights = xi_weights/xi_weights.max()

yw = np.tile(yi_weights, xi_weights.size)
yw = yw.reshape(deproms.T.shape)
xw = np.tile(xi_weights, yi_weights.size)
xw = xw.reshape(deproms.shape)

weights = (xw+yw.T)
weights[weights>1] = 1
weights[bla>100] = 1

sigma = [sigma_y, sigma_x]
weightfilt = sp.ndimage.filters.gaussian_filter(weights, sigma, mode='wrap')
plt.pcolor(weightfilt)

weightsinv = weightfilt *-1 +1

depth1= weightfilt * bla
depth2 = weightsinv * deproms.values

ncroms2 = ncroms.copy()
ncroms2.h.values = depth1+depth2
ncroms2.to_netcdf(romspath + '/sbb_grid_roms_100m.nc')

#
# plt.close('all')
# fig, ax = plt.subplots(ncols=2, nrows=1, sharex=True, sharey=True)
# # ax[0].contourf(lonmerc, latmerc, depmerc, cvalues, cmap=plt.cm.prism,)
# # ax[0].contourf(lonroms, latroms, deproms, cvalues, cmap=plt.cm.prism,)
#
# ax[0].pcolor(lonmerc, latmerc, depmerc, cmap=plt.cm.jet, vmin=2000, vmax=4000)
# ax[0].pcolor(lonroms, latroms, deproms, cmap=plt.cm.jet, vmin=2000, vmax=4000)
#
# ax[0].contour(lonmerc, latmerc, depmerc, cvalues, colors='k',)
# ax[0].contour(lonroms, latroms, deproms, cvalues, colors='w',)
#
#
# # ax[1].contourf(lonmerc, latmerc, depmerc, cvalues, cmap=plt.cm.prism,)
# # ax[1].contourf(lonroms, latroms, depth1 + depth2, cvalues, cmap=plt.cm.prism,)
#
# ax[1].pcolor(lonmerc, latmerc, depmerc, cmap=plt.cm.jet, vmin=2000, vmax=4000)
# ax[1].pcolor(lonroms, latroms, ncroms2.h.values, cmap=plt.cm.jet, vmin=2000, vmax=4000)
#
# ax[1].contour(lonmerc, latmerc, depmerc, cvalues, colors='k',)
# ax[1].contour(lonroms, latroms, ncroms2.h.values, cvalues, colors='w',)
