import xarray as xr
import xarray_decorators
import numpy as np
import xesmf as xe
# from scipy import interpolate as itp
import datetime as dtt
from netCDF4 import date2num
import sys
from utils import utils as ut
import pandas as pd
import gsw
from utils.interpolation import interpolation as itp

class interpolationGlorys2Roms(object):
    def __init__(self, fgrid, fpath, load=False):
        self.dsgrid   = xr.open_dataset(fgrid)
        self.dssource = xr.open_mfdataset(fpath)
        self.dssource  = self.dssource.drop_duplicates(dim='time')

        if load:
            self.dssource.load()

    def monthly_mean(self):
        return self.dssource.groupby('time.month').mean(dim='time')
        

    def set_coords(self, itime, gtype='rho', ztype='rho', time_type=None):
        """
        Set coordinates of glorys to conform with C grid from roms and select
        a given time step of the dataset.  The two options to create the fields
        are monthly means (not climatological) or a given date -- see itime
        argument for more information
        

        itime (integer or datetime):selects the month /date of the interpolation.
        time_type (string):
            'monthly': a monthly mean is calculated and the integer
                  index should be used. For example, if march and may are included in
                  the source netcdf, march and may will be index 0 and 1, respectively.
            None     : uses the date in itime
        gtype (string): ['rho', 'u', 'v']:
            remap horizontal lon, lat to rho, u or v points in Arakawa C-grid
        ztype: ['rho','w']:
            remap vertical positions to rho, or w points in Arakawa C-grid
        """
        dsgrid = self.dsgrid.copy()

        # renaming for interpolation standards ds1.roms.interpolate
        dsgrid = dsgrid.rename({f'lon_{gtype}':'lon', f'lat_{gtype}':'lat'})
        dsgrid = dsgrid.set_coords(['lon', 'lat'])

        # select and set up time
        if time_type is None:
            self.dssource1 = self.dssource
            ds0 = self.dssource1.sel(time=itime, method='nearest')

        elif time_type == 'monthly':
            self.dssource1 = self.monthly_mean()
            self.dssource1 = self.dssource1.rename(month='time')
            ds0 = self.dssource1.isel(time=itime)
        else:
            raise IOError("time_type shoud be None or 'monthly'")
        
        self.ds1 = ds0.roms.interpolate(dsgrid,
                {'longitude':'lon', 'latitude':'lat'},
                xesmf_regridder_kwargs=dict(extrap_method='nearest_s2d'))
        a = self.ds1

        # glorys interpolated horizontally to roms
        self.zz   = -np.tile(a.depth.values[:,None,None], (1, *a.lon.shape) )
        self.lonz = np.broadcast_to(a.lon.values, (ds0.depth.size, *a.lon.shape))
        self.latz = np.broadcast_to(a.lat.values, (ds0.depth.size, *a.lat.shape))


        self.zrho   = self.dsgrid.roms.s2z(ztype).transpose(f's_{ztype}',f'eta_{ztype}', f'xi_{ztype}')
        
        self.zu     =  self.zrho#.rename(lon_rho='lon_u', lat_rho='lat_u', eta_rho='eta_u', xi_rho='xi_u')
        self.zu     = self.zu[:,:,:-1]
        # self.zu.coords['lon_u'].values = self.dsgrid.lon_u
        # self.zu.coords['lat_u'].values = self.dsgrid.lat_u

        self.zv     =  self.zrho#.rename(lon_rho='lon_v', lat_rho='lat_v', eta_rho='eta_v', xi_rho='xi_v')
        self.zv     = self.zv[:,:-1,:]
        # self.zv.coords['lon_v'].values = self.dsgrid.lon_v
        # self.zv.coords['lat_v'].values = self.dsgrid.lat_v

        
        if gtype=='u':
            self.zg = self.zu
            self.long = np.broadcast_to(dsgrid['lon'].values, self.zu.shape)
            self.latg = np.broadcast_to(dsgrid['lat'].values, self.zu.shape)

        elif gtype=='v':
            self.zg = self.zv
            self.long = np.broadcast_to(dsgrid['lon'].values, self.zv.shape)
            self.latg = np.broadcast_to(dsgrid['lat'].values, self.zv.shape)

        elif gtype=='rho':
            self.zg = self.zrho
            self.long = np.broadcast_to(dsgrid['lon'].values, self.zrho.shape)
            self.latg = np.broadcast_to(dsgrid['lat'].values, self.zrho.shape)

        else:
            IOError('gtype must be u,v or rho')


        return self.dssource1['time']


    def calculate_geostrophy(self):
        ds = self.dssource
        xm, ym = np.meshgrid(ds.latitude.values,
                            ds.longitude.values,)

        f = gsw.f(ym)
        g = 9.8



        ds['f'] = (['latitude', 'longitude'], f)
        ds['psi'] = ds['zos'] * 9.8/ ds['f']

        ds['dpsi_x'] = (['time', 'latitude', 'longitude'], np.gradient(ds['psi'].values, axis=2))
        ds['dpsi_y'] = (['time', 'latitude', 'longitude'], np.gradient(ds['psi'].values, axis=1))
        ds.drop('psi')


        # calcualte distances
        lons = ds.longitude.values
        lats = ds.latitude.values

        # provisory dy dx
        dy = gsw.distance(xm,ym,axis=1)
        dx = gsw.distance(xm,ym,axis=0)

        # cumulative sum 
        x = np.cumsum(dx, axis=1)
        y = np.cumsum(dy, axis=0)

        # recalculating distances to preserve dimensions
        x = np.insert(x, 0, values=0, axis=0)
        y = np.insert(y, 0, values=0, axis=1)

        dx = np.gradient(x, axis=1)
        dy = np.gradient(y, axis=0)

        # set in xarray dataset 
        ds['dx'] = (['latitude', 'longitude'], dx)
        ds['dy'] = (['latitude', 'longitude'], dy)

        # calculate geostrophic velocities
        ds['ug'] = -ds['dpsi_y']/ds['dy']
        ds['vg'] = ds['dpsi_x']/ds['dx']


        # plt.close('all')
        # fig = plt.figure(figsize=[20,20])
        # plt.contourf(xm,ym,ds['zos'][0], levels=30,zorder=0, cmap=plt.cm.jet)
        # plt.quiver(xm[::2,::2],ym[::2,::2], ds['ug'][0,::2,::2], ds['vg'][0,::2,::2], scale=20, width=0.0005)
        # plt.savefig('bla.png', dpi=300)

    def interpolate3d(self,varb):
        zs = self.zg
        outfield = np.zeros(self.zg.shape)

        xm = self.lonz
        ym = self.latz
        zm = self.zz[::-1]
        varb = self.ds1[varb]
        varb = varb.bfill('xi_rho')
        varb = varb.ffill('xi_rho')
        varb = varb.bfill('eta_rho')
        varb = varb.ffill('eta_rho')

        varb = varb.values[::-1]

        nzs = zs.shape[0]
        nz, ny, nx = xm.shape
        aux = itp.interp3d_along_axis0(zs,
                                    xm,ym,zm,
                                    varb,
                                    nzs,nx,ny,nz)
        aux = np.where(aux==9999.0, np.nan, aux)
        outfield[...] = aux
        return outfield

    # def interpolate3d(self, varb):
    #     outfield = np.zeros(self.zg.shape)
    #     for i in range(self.zg.shape[2]):
    #         print(i)
    #         zsurfi   = self.zg[:,:,i].values
    #         lonsurfi = self.long[:,:,i]
    #         values = self.dsglo[varb][:,:,i].values.ravel()

    #         zsurf   = self.zz[:,:,i]
    #         lonsurf = self.lonz[:,:,i]


    #         aux = itp.griddata((lonsurf.ravel(),zsurf.ravel()),
    #                         values,
    #                         (lonsurfi.ravel(), zsurfi.ravel())
    #                             )
    #         aux = aux.reshape(zsurfi.shape)

    #         outfield[:,:,i] = aux
    #     return outfield


def nearest_interpolation(dsgrid, varb, hgrid='rho'):
    dsgrid[varb] = dsgrid[varb].bfill('s_rho')
    dsgrid[varb] = dsgrid[varb].ffill('s_rho')
    dsgrid[varb] = dsgrid[varb].ffill(f'xi_{hgrid}')
    dsgrid[varb] = dsgrid[varb].bfill(f'xi_{hgrid}')
    dsgrid[varb] = dsgrid[varb].ffill(f'eta_{hgrid}')
    dsgrid[varb] = dsgrid[varb].bfill(f'eta_{hgrid}')
    return dsgrid


if __name__ == '__main__':


    if len(sys.argv) > 1:
        reference = sys.argv[1]
    else:
        reference = 'ceresIV_2.012'

    dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
    dicts = dicts[reference]
    grid  = dicts['grid_dir']
    source= dicts['clim.src_file']
    output= dicts['clim.outfile']
    tstart = dicts['clim.date'][0]
    tfinal = dicts['clim.date'][1]
    tfreq  = dicts['clim.date'][2]
    daterange = pd.date_range(start=tstart, end=tfinal, freq=tfreq)

    print('starting loop')
    for datetime in daterange:
        print(datetime)
        dsgrid = xr.open_dataset(grid)
        tref = datetime.to_pydatetime()
        tref1 = date2num(tref, 'days since 1990-01-01 00:00:00')


        # glorys data in interpolation Glorys2Roms is read with xr.open_mfdataset. 
        # If using load kwarg, the script will run faster but RAM will be used more intesively
        # open_mfdataset may be slow depending on the size of the grid
        q = interpolationGlorys2Roms(grid,
                                    source,
                                    load=True)
        q.calculate_geostrophy()
        q.set_coords(tref, gtype='rho', time_type=None)
        temp = q.interpolate3d('thetao')
        salt = q.interpolate3d('so')
        # zeta = q.interpolate3d('zos')

        # interpolat3d u and v onto rho grid (due to rotation)
        q.set_coords(tref, gtype='rho', time_type=None)
        v  = q.interpolate3d('vo')
        t = q.set_coords(tref, gtype='rho', time_type=None)
        u    = q.interpolate3d('uo')

        print('interpolation done')

        # rotating u and v as eta and xi coordinate components
        rot = dsgrid.angle.values

        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v,u)

        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)



        # interpolate u and v onto rho grid (due to rotation)
        v  = q.ds1['ug']
        u  = q.ds1['vg']

        # rotating u and v as eta and xi coordinate components
        rot = dsgrid.angle.values

        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v,u)

        ug = -mag * np.sin(angle + rot)
        vg = mag * np.cos(angle + rot)

        dsgrid = dsgrid.assign_coords(time=[tref])
        dsgrid = dsgrid.assign_coords(temp_time=[tref])
        dsgrid = dsgrid.assign_coords(salt_time=[tref])

        dsgrid['temp'] = (('temp_time', 's_rho', 'eta_rho', 'xi_rho'), [temp])
        dsgrid['salt'] = (('salt_time', 's_rho', 'eta_rho', 'xi_rho'), [salt])
        dsgrid['u'] = (('time', 's_rho', 'eta_u', 'xi_u'), [u1[:,:,:-1]])
        dsgrid['v'] = (('time', 's_rho', 'eta_v', 'xi_v'), [v1[:,:-1,:]])

        dsgrid['ubar'] = (('time', 'eta_u', 'xi_u'), [ug[:,:-1]])
        dsgrid['vbar'] = (('time', 'eta_v', 'xi_v'), [vg[:-1,:]])


        for v,c in zip(['temp', 'salt','v','u'],['rho','rho','v','u']):
            dsgrid = nearest_interpolation(dsgrid, v, hgrid=c)


        dsgrid['time'].attrs = {}
        dsgrid['time'].attrs['long_name'] = 'time'
        # dsgrid['time'].attrs['unit'] = 'day'
        # dsgrid['time'].attrs['cycle_length'] = np.array(12).astype(float)

        dsgrid.to_netcdf(output % str(tref).replace(' ','T'))
        print(output % str(tref).replace(' ','T') + ' saved')
print('Done')


