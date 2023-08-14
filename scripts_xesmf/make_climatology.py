import xarray as xr
import xarray_decorators
import numpy as np
import xesmf as xe
from scipy import interpolate as itp
import datetime as dtt
from netCDF4 import date2num
import sys
from utils import utils as ut

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

        if time_type is None:
            self.dssource1 = self.dssource
            ds0 = self.dssource1.sel(time=itime, method='nearest')

        elif time_type == 'monthly':
            self.dssource1 = self.monthly_mean()
            self.dssource1 = self.dssource1.rename(month='time')
            ds0 = self.dssource1.isel(time=itime)
        else:
            raise IOError("time_type shoud be None or 'monthly'")
        
        ds1 = ds0.roms.interpolate(dsgrid, {'longitude':'lon', 'latitude':'lat'})
        self.ds1= ds1.roms.interpolate(dsgrid, gtype)

        a = ds1.roms.interpolate(dsgrid, 'bla')
        self.dsglo = a


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

    def interpolate(self, varb):
        outfield = np.zeros(self.zg.shape)
        for i in range(self.zg.shape[2]):
            print(i)
            zsurfi   = self.zg[:,:,i].values
            lonsurfi = self.long[:,:,i]
            values = self.dsglo[varb][:,:,i].values.ravel()

            zsurf   = self.zz[:,:,i]
            lonsurf = self.lonz[:,:,i]


            aux = itp.griddata((lonsurf.ravel(),zsurf.ravel()),
                            values,
                            (lonsurfi.ravel(), zsurfi.ravel())
                                )
            aux = aux.reshape(zsurfi.shape)

            outfield[:,:,i] = aux
        return outfield


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
        reference = 'ceresIV_2.9'

    dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
    dicts = dicts[reference]
    grid  = dicts['grid_dir']
    source= dicts['clim.src_file']
    output= dicts['clim.outfile']



    for datetime in dicts['clim.date']:
        dsgrid = xr.open_dataset(grid)
        tref = dtt.datetime(*datetime)
        tref1 = date2num(tref, 'days since 1990-01-01 00:00:00')


        # glorys data in interpolation Glorys2Roms is read with xr.open_mfdataset. 
        # If using load kwarg, the script will run faster but RAM will be used more intesively
        # open_mfdataset may be slow depending on the size of the grid
        q = interpolationGlorys2Roms(grid,
                                    source,
                                    load=True)


        q.set_coords(tref, gtype='rho', time_type=None)
        temp = q.interpolate('thetao')
        salt = q.interpolate('so')
        # zeta = q.interpolate('zos')

        # interpolate u and v onto rho grid (due to rotation)
        q.set_coords(tref, gtype='rho', time_type=None)
        v  = q.interpolate('vo')
        t = q.set_coords(tref, gtype='rho', time_type=None)
        u    = q.interpolate('uo')

        # rotating u and v as eta and xi coordinate components
        rot = dsgrid.angle.values

        mag = (u**2 + v**2)**0.5
        angle = np.arctan2(v,u)

        u1 = -mag * np.sin(angle + rot)
        v1 = mag * np.cos(angle + rot)


        dsgrid = dsgrid.assign_coords(time=[tref])
        dsgrid = dsgrid.assign_coords(temp_time=[tref])
        dsgrid = dsgrid.assign_coords(salt_time=[tref])

        dsgrid['temp'] = (('temp_time', 's_rho', 'eta_rho', 'xi_rho'), [temp])
        dsgrid['salt'] = (('salt_time', 's_rho', 'eta_rho', 'xi_rho'), [salt])
        dsgrid['u'] = (('time', 's_rho', 'eta_u', 'xi_u'), [u1[:,:,:-1]])
        dsgrid['v'] = (('time', 's_rho', 'eta_v', 'xi_v'), [v1[:,:-1,:]])


        for v,c in zip(['temp', 'salt','v','u'],['rho','rho','v','u']):
            dsgrid = nearest_interpolation(dsgrid, v, hgrid=c)




        dsgrid['time'].attrs = {}
        dsgrid['time'].attrs['long_name'] = 'time'
        # dsgrid['time'].attrs['unit'] = 'day'
        # dsgrid['time'].attrs['cycle_length'] = np.array(12).astype(float)

        dsgrid.to_netcdf(output % str(tref).replace(' ','T'))

