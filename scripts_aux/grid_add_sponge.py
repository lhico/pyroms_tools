import xarray as xr
import numpy as np

nc = xr.open_dataset('bacia_santos.nc')
nc['visc_factor'] = (['eta_rho', 'xi_rho'], np.zeros(nc.mask_rho.shape))
nc['visc_factor'].attrs['long_name'] = "horizontal viscosity sponge factor" ;
nc['visc_factor'].attrs['valid_min'] = 0. ;
nc['visc_factor'].attrs['coordinates'] = "lon_rho lat_rho" ;


nc['diff_factor'] = (['eta_rho', 'xi_rho'], np.zeros(nc.mask_rho.shape))
nc['diff_factor'].attrs['long_name'] = "horizontal diffusivity sponge factor" ;
nc['diff_factor'].attrs['valid_min'] = 0. ;
nc['diff_factor'].attrs['coordinates'] = "lon_rho lat_rho" ;

for i in ['visc_factor', 'diff_factor']:
    # for cont in range(30,-2,-1):
    #     nc[i].values[cont-1:-cont+1,-cont] =  4
    #     nc[i].values[cont-1,cont-1:-cont+1] =  4
    #     nc[i].values[-cont,cont-1:-cont+1] =  4


    nc[i].values[-5:,:] = 1e1
    nc[i].values[:,-5:] = 1e1
    nc[i].values[:5,:] = 1e1

nc.to_netcdf('bacia_santos1.nc')