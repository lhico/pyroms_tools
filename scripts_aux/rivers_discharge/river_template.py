from utils import maskedge
from utils import click_map as cp
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from utils import utils as ut
import pandas as pd
import sys


def scatter(event, tracker):
    cont = tracker.cont -1
    print(tracker.points[cont])
    j = tracker.points[cont]['j']
    i = tracker.points[cont]['i']

    tracker.ax.scatter(tracker.dataset['lon_rho'][j,i],
                       tracker.dataset['lat_rho'][j,i],
                       transform=tracker.proj,
                       c='b', marker='x')
    fig.canvas.draw_idle()


def river_points(points, start='2020-01-01', end='2020-05-01', s_rho=40):
    points = pd.DataFrame(points)
    points = points.T
    print(points)


    value = 10
    nzlevel = s_rho.size
    nrivers = points.shape[0]
    print(points.shape)
    t = pd.date_range(start=start, end=end)
    aux = np.ones(t.size)

    riverid = np.arange(nrivers)
    time      = t.copy()
    xpos      = points['i'].values
    epos      = points['j'].values
    lon      = points['lon'].values
    lat      = points['lat'].values
    sign      = points['sign'].values
    direction = points['direc'].values
    Vshape    = np.tile(np.ones(nzlevel)/nzlevel, (nrivers  ,1)).T
    transport = np.tile(aux * value, (nrivers,1)).T
    temp      = np.ones([t.size, nzlevel, nrivers]) * 24
    salt      = np.ones([t.size, nzlevel, nrivers]) * 0
    print(Vshape)
    ncout = ut.river_netcdf(riverid, s_rho, time,
                        xpos, epos, sign, direction,
                        Vshape, transport, temp, salt,
                        lon, lat)
    return ncout


if __name__ == '__main__':
    """
    Create a river file used by ROMS.
    First a map is outputted, where the edge points are marked by red
    dots. On the map, we can select the dots as needed (a blue cross
    will be included in selected points)

    After selecting a point we have the option to include the direction
    in the roms staggered grid by pressing the keyboard keys 0, 1, or 2
    (related to u, v, rho points), and + (0) for downstream (upstream)
    currents relative to u.

    when the selection is finished, press q.  At this point, the netcdf
    template should be ready for use.
    """
    ffile = sys.argv[1] if len(sys.argv) == 4 else "modified_grid.nc"
    start=sys.argv[2] if len(sys.argv) == 4 else '2010-01-01'
    end=sys.argv[3] if len(sys.argv) == 4 else '2023-01-01'

    points = maskedge.main(ffile)
    points = np.array(points)
    nc = xr.open_dataset(ffile)

    # 0 is u, 1 is v
    # +1 (-1) diverges (converges) from left to right in u direction
    plt.ion()
    plt.close('all')
    tracker = cp.IndexTracker(nc, ttype='get_index')
    while tracker.wait:
        # tracker = cp.IndexTracker(nc, ttype='mask')
        fig, ax = tracker.get_figure()
        x = tracker.dataset["lon_rho"].values[points[:,1], points[:,0]]
        y = tracker.dataset["lat_rho"].values[points[:,1], points[:,0]]
        ax.scatter(x, y, transform=tracker.proj, zorder=10, c='r', s=1)

        tracker.ax.figure.canvas.mpl_connect('button_press_event', lambda b: ax.scatter(x, y, zorder=10, c='r', s=1))
        tracker.ax.figure.canvas.mpl_connect('button_press_event', lambda b: scatter(b, tracker))

        plt.pause(2)

        fig.canvas.draw_idle()
        # fig.tight_layout()
        tracker.show()
        plt.pause(2)

    ncout = river_points(tracker.points, start=start, end=end, s_rho=nc.s_rho.values)
    ncout.to_netcdf('river_template.nc')
