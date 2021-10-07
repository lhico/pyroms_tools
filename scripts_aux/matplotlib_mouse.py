import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def onclick(event):
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

    X[i,j] = 0 if X[i,j] ==1 else 1

    # fig.imshow(X)
    ax.clear()
    plt.gcf().canvas.draw_idle()
    plot(maskstr=masklist[ivar], lonstr=lonlist[ivar], latstr=latlist[ivar])


def nearest_index(ds, lat, lon, latstr='lat_rho', lonstr='lon_rho'):
    # Find absolute difference between requested point and the grid coordinates.
    abslat = np.abs(ds[latstr] - lat)
    abslon = np.abs(ds[lonstr] - lon)

    # Create grid of the maximum values of the two absolute grids
    c = (abslon**2 + abslat**2)**0.5

    # Find location where lat/lon minimum absolute value intersects
    y, x = np.where(c == np.min(c))
    return y,x


def plot(maskstr='mask_rho', lonstr='lon_rho', latstr='lat_rho'):
    ncchange[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.1, alpha=0.5, cmap=plt.cm.jet)
    ncref[maskstr].plot(x=lonstr, y=latstr, add_colorbar=False, edgecolors='k', linewidths=0.1, alpha=0.5, cmap=plt.cm.jet)
    plt.axis('equal')
    plt.xlim(np.min(ncchange.lon_rho.values), np.max(ncchange.lon_rho.values))
    plt.ylim(np.min(ncchange.lat_rho.values), np.max(ncchange.lat_rho.values))


masklist = ['mask_rho', 'mask_u', 'mask_v']
lonlist  = ['lon_rho', 'lon_u', 'lon_v']
latlist  = ['lat_rho', 'lat_u', 'lat_v']

ncref = xr.open_dataset('bacia_santos_donord.nc')
ncchange = xr.open_dataset('bacia_santos_ref1d.nc')


for ivar in range(3):
    
    ncref, ncchange = map(lambda x: x.set_coords([lonlist[ivar], latlist[ivar]]), [ncref, ncchange])

    fig, ax = plt.subplots()
    X = ncchange[masklist[ivar]].values
    Y = ncref[masklist[ivar]].values
    plot(maskstr=masklist[ivar], lonstr=lonlist[ivar], latstr=latlist[ivar])


    plt.gcf().canvas.mpl_connect('button_press_event', onclick)

    plt.show()

ncchange.to_netcdf('test.nc')