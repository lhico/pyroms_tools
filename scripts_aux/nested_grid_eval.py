import xarray as xr
import matplotlib.pyplot as plt
from dry_package.plot_schemes import maps
import numpy as np
import glob
import cartopy.feature as cfeature

land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor='0.5')

def u_eastward(tp='pcolor'):
    if tp == 'contourf':
        ax = maps.make_overall()#plt.figure()
        cf =nc1.u_eastward[-1,-1].plot.contourf(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, levels=np.arange(-0.2,0.21,0.01), cmap=plt.cm.RdBu_r, zorder=1)
        nc2.u_eastward[-1,-1,1:-1,1:-1].plot.contourf(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, levels=np.arange(-0.2,0.21,0.01), cmap=plt.cm.RdBu_r, zorder=3)
        nc3.u_eastward[-1,-1,1:-1,1:-1].plot.contourf(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, levels=np.arange(-0.2,0.21,0.01), cmap=plt.cm.RdBu_r, zorder=5)
    else:
        ax = maps.make_overall()#plt.figure()
        cf = nc1.u_eastward[-1,-1].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.5, vmax=0.5, cmap=plt.cm.RdBu_r, zorder=1)
        nc2.u_eastward[-1,-1,1:-1,1:-1].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.5, vmax=0.5, cmap=plt.cm.RdBu_r, zorder=3)
        nc3.u_eastward[-1,-1,1:-1,1:-1].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.5, vmax=0.5, cmap=plt.cm.RdBu_r, zorder=5)

    nc1.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(0,100,2.5), cmap=plt.cm.jet, linewidths=0.5, zorder=2)
    nc2.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(0,100,2.5), cmap=plt.cm.jet, linewidths=0.5, zorder=4, linestyles='--')
    nc3.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(0,100,2.5), cmap=plt.cm.jet, linewidths=0.5, zorder=6, linestyles=':')
    
    ax.add_feature(land_10m, zorder=10)
    plt.colorbar(cf)


def mask_rho():
    ax = maps.make_overall()#plt.figure()
    nc1.mask_rho.plot(x='lon_rho', y='lat_rho', add_colorbar=False, edgecolors='k', alpha=0.5, linewidth=0.1, ax=ax)
    nc2.mask_rho[2:-2,2:-2].plot(x='lon_rho', y='lat_rho', add_colorbar=False, edgecolors='k', alpha=0.5, linewidth=0.1, ax=ax)
    nc3.mask_rho.plot(x='lon_rho', y='lat_rho', add_colorbar=False, edgecolors='k', alpha=0.5, linewidth=0.1, ax=ax)


def h():
    ax = maps.make_overall()#plt.figure()
    kwargs = dict(add_colorbar=False, edgecolors='k', alpha=1,
                  linewidth=0.05, ax=ax, cmap=plt.cm.nipy_spectral,
                  vmin=5, vmax=80)
    nc1.h.plot(x='lon_rho', y='lat_rho', **kwargs)
    nc2.h.plot(x='lon_rho', y='lat_rho', **kwargs)
    nc3.h.plot(x='lon_rho', y='lat_rho', **kwargs)

def salt():
    ax = maps.make_overall()#plt.figure()
    nc1.salt[-1,0].plot(x='lon_rho', y='lat_rho', vmin=34, vmax=36.6, cmap=plt.cm.jet, add_colorbar=False, ax=ax)
    nc2.salt[-1,0].plot(x='lon_rho', y='lat_rho', vmin=34, vmax=36.6, cmap=plt.cm.jet, add_colorbar=False, ax=ax)
    nc3.salt[-1,0].plot(x='lon_rho', y='lat_rho', vmin=34, vmax=36.6, cmap=plt.cm.jet, add_colorbar=False, ax=ax)


if __name__ == '__main__':
    # nc = ['bacia_santos_donor.nc', 'nestref1_glorys_ic2.nc', 'bacia_santos_ref2b.nc']
    # nc = ['input2/bacia_santos_donorc.nc', 'input2/bacia_santos_ref1b.nc', 'input2/bacia_santos_ref2b.nc']
    nc = ['roms_avg.nc', 'roms_avg_nest.nc', 'roms_avg_nest2.nc']
    nc = map(xr.open_dataset, nc)
    nc1,nc2,nc3 = map(lambda x: x.set_coords(['lon_rho','lat_rho']), nc)

    # plt.figure()
    # nc1.h.plot(x='lon_rho', y='lat_rho', vmax=80, cmap=plt.cm.jet, add_colorbar=False)
    # nc2.h.plot(x='lon_rho', y='lat_rho', vmax=80, cmap=plt.cm.jet, add_colorbar=False)
    # nc3.h.plot(x='lon_rho', y='lat_rho', vmax=80, cmap=plt.cm.jet, add_colorbar=False)


    # plt.figure()
    # nc1.zeta[-1].plot(x='lon_rho', y='lat_rho', vmin=-0.5, vmax=0.5, cmap=plt.cm.jet, add_colorbar=False)
    # nc2.zeta[-1].plot(x='lon_rho', y='lat_rho', vmin=-0.5, vmax=0.5, cmap=plt.cm.jet, add_colorbar=False)
    # nc3.zeta[-1].plot(x='lon_rho', y='lat_rho', vmin=-0.5, vmax=0.5, cmap=plt.cm.jet, add_colorbar=False)


    # ax = maps.make_overall()#plt.figure()
    # nc1.u_eastward[-1,-1].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.2, vmax=0.2, cmap=plt.cm.RdBu_r)
    # nc2.u_eastward[-1,-1,2:-2,2:-2].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.2, vmax=0.2, cmap=plt.cm.RdBu_r)
    # nc3.u_eastward[-1,-1,2:-2,2:-2].plot(x='lon_rho', y='lat_rho', add_colorbar=False, linewidth=0.1, vmin=-0.2, vmax=0.2, cmap=plt.cm.RdBu_r)


