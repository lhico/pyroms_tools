import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from dry_package.grid import interp
from dry_package.decorators import decorate as dc
from matplotlib.backend_bases import MouseButton

import threading
import queue

class IndexTracker:
    def __init__(self, dataset, ttype='mask'):
        self.dataset = dataset
        self.points  = {}
        self.cont    = 0
        self.wait    = True

        # Create a cartopy projection
        self.proj = ccrs.PlateCarree()
        x= self.proj.transform_points(self.proj, dataset.lon_rho.values, dataset.lat_rho.values, )
        y= x[...,1]
        x= x[...,0]
        self.dataset[['mask_rho', 'mask_u', 'mask_v', 'mask_psi']].load()

        self.dataset['lon_rho'].values = x 
        self.dataset['lat_rho'].values = y
        self.curvgrid = interp.curvilinear_grid_op(x,y)
        self.fig, self.ax = plt.subplots(subplot_kw={"projection": self.proj})
        self.ax.set_title('Click on a point')

        # Use pcolormesh with cartopy
        self.mesh = self.ax.pcolormesh(
            self.dataset.lon_rho, self.dataset.lat_rho, self.dataset["mask_rho"],
            transform=self.proj, shading='auto', alpha=0.5
        )
        self.ax.set_extent([
            np.min(self.dataset.lon_rho), np.max(self.dataset.lon_rho),
            np.min(self.dataset.lat_rho), np.max(self.dataset.lat_rho)
        ])

        self.ax.coastlines(resolution='10m')
        self.coords = []

        # keyboard and buttom inputs with options whe needed
        # the options are in on_click and on_press methods
        self.ax.figure.canvas.mpl_connect('button_press_event', lambda x: self.onclick(x, ttype=ttype))
        self.ax.figure.canvas.mpl_connect('key_press_event', lambda x: self.onpress(x, ttype))

    def get_figure(self):
        return self.fig, self.ax

    @dc.log_time
    def onclick(self, event, ttype='mask'):
        if event.button is MouseButton.RIGHT:
            if ttype == 'mask':
                return self._onclick_mask(event)
            elif ttype == 'get_index':
                return self._onclick_get_index(event)
            self.ax.figure.canvas.mpl_disconnect(self.cd1)

    def _onclick_get_index(self, event):
        if event.inaxes == self.ax:
            x, y = event.xdata, event.ydata
            i, j = self.get_index(x,y)
            self.points[self.cont] = {'i': None, 'j':None, 'lat':None, 'lon':None}
            self.points[self.cont]['j'] = j
            self.points[self.cont]['i'] = i
            self.points[self.cont]['lat'] = y
            self.points[self.cont]['lon'] = x
            print(f'(j,i)={j,i}, (y,x)={y,x}')
            self.cont += 1

    def _onclick_mask(self, event):
        if event.inaxes == self.ax:
            x, y = event.xdata, event.ydata
            i, j = self.get_index(x,y)
            # i, j = self.get_index(x, y)
            mask_rho = self.dataset['mask_rho']

            print(mask_rho.values[j,i])
            if mask_rho.values[j,i]:
                mask_rho.values[j,i] = 0
                self.dataset.h[j,i] = 5.
                self.dataset.hraw[:,j,i] = 5.
            else:
                mask_rho.values[j,i] = 1
                self.dataset.h.values[j,i] = 10.
                self.dataset.hraw[:,j,i] = 10.
            # self.ax.clear()
            self.mesh = self.ax.pcolormesh(
                self.dataset.lon_rho, self.dataset.lat_rho, self.dataset["mask_rho"],
                transform=self.proj, alpha=0.5
            )
            self.ax.coastlines(resolution='10m')
            self.fig.canvas.draw_idle()

            # self.dataset['mask_u'].values = mask_rho[:,:-1].values * mask_rho[:,1:].values
            # self.dataset['mask_v'].values = mask_rho[:-1,:].values * mask_rho[1:,:].values
            # self.dataset['mask_psi'].values = self.dataset.mask_u[:-1, :].values * \
            #                                   self.dataset.mask_v[:, :-1].values

            # self.dataset['mask_psi'][j,i] = mask_rho[j,i]




    def onpress(self, event, ttype='mask'):
        if ttype == 'mask':
            return self._onpress_mask(event)
        elif ttype == 'get_index':
            return self._onpress_get_index(event)
        self.ax.figure.canvas.mpl_disconnect(self.cd2)

    def _onpress_get_index(self, event):
        if event.key == "z":
            if self.points:
                remove = self.points.pop(self.cont - 1)
                print(f"removed {remove}")
                self.cont -= 1
            elif self.cont == 0:
                print("No clicks to undo.")
        elif event.key == "q":
            self._onpress_close()
            self.wait = False
        elif event.key in ['0','1','2']:
            print(event.key)
            self.points[self.cont-1]['direc'] = event.key
            print(f'{self.points[self.cont-1]}')
        elif event.key in ['+','-']:
            print(event.key)
            sign = {'+': 1, '-': -1}
            self.points[self.cont-1]['sign'] = sign[event.key]
            print(f'{self.points[self.cont-1]}')
        elif event.key == 'c':
            self.ax.clear()
            self.fig.canvas.draw_idle()



    def _onpress_mask(self, event):

        if event.key == "z":
            if self.coords:
                self.coords.pop()
                print("Last click undone.")
            else:
                print("No clicks to undo.")
        elif event.key == "q":
            self._onpress_close()
            self.wait = False
        elif event.key == 'u':
            fig,ax = plt.subplots()
            ax.pcolor(self.dataset['lon_u'], self.dataset['lat_u'], self.dataset['mask_u'])
        elif event.key == 'v':
            fig,ax = plt.subplots()
            ax.pcolor(self.dataset['lon_v'], self.dataset['lat_v'], self.dataset['mask_v'])
        # elif event.key == 'g':
        #     self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

  

    def user_input(self):
        return input("Enter something: ")

    def _onpress_close(self):
        self.fig.clf()


    def on_motion(self, event):
        if event.inaxes is not None:
            x, y = event.xdata, event.ydata
            i, j = self.get_index(x,y)
            self.ax.format_coord = (
                lambda x, y: f"Lon: {x:.2f}  Lat: {y:.2f}  Index: ({i}, {j})"
            )

    def get_index(self, x, y):
        _, ji, _ = self.curvgrid.query(x,y,k=2)
        # lat_index = (
        #     self.dataset["lat_rho"]
        #     .where(
        #         (self.dataset["lat_rho"] <= y) & (self.dataset["lon_rho"] <= x),
        #         drop=True,
        #     )
        #     .argmax(dim="eta_rho")
        # )
        # lon_index = (
        #     self.dataset["lon_rho"]
        #     .where(
        #         (self.dataset["lat_rho"] <= y) & (self.dataset["lon_rho"] <= x),
        #         drop=True,
        #     )
        #     .argmax(dim="xi_rho")
        # )
        # j = list(set(lat_index.values))[0]
        # i = list(set(lon_index.values))[0]

        return ji[0][1], ji[0][0]

    def show(self):
        self.fig.show()


if __name__ == '__main__':
    # Example usage
    nc = xr.open_dataset("sbb_grid_roms.nc")

    plt.close('all')
    # tracker = IndexTracker(nc, ttype='get_index')
    tracker = IndexTracker(nc, ttype='mask')
    tracker.show()
