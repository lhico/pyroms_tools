import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.backend_bases import MouseButton

def onclick(event):
    global i, j

    print('Returned indices')
    print(event.xdata, event.ydata)
    print('mapping back:')
    x = event.xdata
    y = event.ydata

    tx = f"lat: {y}, lon: {x}"
    print(tx)

    i, j = nearest_index(ncchange, y, x)
    tx = f"i: {i}, j: {j}, Value: {X[i,j]}"
    print(tx)
    
    # fig.imshow(X)
    ax.clear()
    plt.gcf().canvas.draw_idle()
    plot(maskstr=masklist[ivar], lonstr=lonlist[ivar], latstr=latlist[ivar], ref='donor')
    

    if event.button is MouseButton.RIGHT:
        print('disconnecting callback')
        plt.gcf().canvas.mpl_disconnect(cid)

def nearest_index(ds, lat, lon, latstr='lat_rho', lonstr='lon_rho'):
    # Find absolute difference between requested point and the grid coordinates.
    abslat = np.abs(ds[latstr] - lat)
    abslon = np.abs(ds[lonstr] - lon)

    # Create grid of the maximum values of the two absolute grids
    c = (abslon**2 + abslat**2)**0.5

    # Find location where lat/lon minimum absolute value intersects
    y, x = np.where(c == np.min(c))
    return y,x


def plot(maskstr='mask_rho', lonstr='lon_rho', latstr='lat_rho', vmin=5, vmax=15, ref='donor'):
    if ref=='donor':
        cb = ncref[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.3, alpha=0.5, cmap=plt.cm.nipy_spectral, vmin=vmin, vmax=vmax)
        ncchange[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.3, alpha=0.5, cmap=plt.cm.nipy_spectral, vmin=vmin, vmax=vmax)
        plt.axis('equal')

        plt.xlim(np.min(ncchange.lon_rho.values), np.max(ncchange.lon_rho.values)-0.2)
        plt.ylim(np.min(ncchange.lat_rho.values), np.max(ncchange.lat_rho.values)-0.15)
        # plt.xlim(-45.2,-44.8)
        # plt.ylim(-23.8,-23.4)

        plt.xlim(-45.65,-45.55)
        plt.ylim(-23.77,-23.86)

    else:
        ncchange[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.1, alpha=0.5, cmap=plt.cm.nipy_spectral, vmin=vmin, vmax=vmax)
        cb = ncref[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.1, alpha=0.5, cmap=plt.cm.nipy_spectral, vmin=vmin, vmax=vmax)
        plt.axis('equal')
        plt.xlim(np.min(ncref.lon_rho.values), np.max(ncref.lon_rho.values))
        plt.ylim(np.min(ncref.lat_rho.values), np.max(ncref.lat_rho.values)-0.15)
        # plt.xlim(-45.2,-44.8)
        # plt.ylim(-23.8,-23.4)

        plt.xlim(-45.65,-45.55)
        plt.ylim(-23.77,-23.86)

    cl1 = ncchange.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(vmin,vmax,1), cmap=plt.cm.jet)
    cl2 = ncref.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(vmin,vmax,1), cmap=plt.cm.jet, linestyles='--')
    plt.clabel(cl1, colors='k', fontsize=5)
    plt.clabel(cl2, colors='k', fontsize=5)


        

    # plt.colorbar(cb)

def on_move(event):
    # get the x and y pixel coords
    x, y = event.x, event.y


masklist = ['h', 'mask_u', 'mask_v']
lonlist  = ['lon_rho', 'lon_u', 'lon_v']
latlist  = ['lat_rho', 'lat_u', 'lat_v']

ncchange = xr.open_dataset('bacia_santos_ref2c.nc')
ncref = xr.open_dataset('bacia_santos_ref1d.nc')

plt.close('all')
q = None
for ivar in range(1):
    
    ncref, ncchange = map(lambda x: x.set_coords([lonlist[ivar], latlist[ivar]]), [ncref, ncchange])

    fig, ax = plt.subplots()
    X = ncchange[masklist[ivar]].values
    Y = ncref[masklist[ivar]].values
    plot(maskstr=masklist[ivar], lonstr=lonlist[ivar], latstr=latlist[ivar], ref='donor')

    while q != 'out':
        plt.ion()
        plt.show()
        cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
        plt.ioff()

        q = input("input value (to quit write 'out'): ")
        X[i,j] = float(q)
        ncchange[masklist[ivar]].values = X
    

ncchange.to_netcdf('test.nc')