import xarray as xr
import numpy as np

nc = xr.open_dataset('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_grid01_2.012_smooth.nc')


nc['visc_factor'] = (['eta_rho', 'xi_rho'], np.ones(nc.mask_rho.shape))
nc['visc_factor'].attrs['long_name'] = "horizontal viscosity sponge factor" ;
nc['visc_factor'].attrs['valid_min'] = 0. ;
nc['visc_factor'].attrs['coordinates'] = "lon_rho lat_rho" ;


nc['diff_factor'] = (['eta_rho', 'xi_rho'], np.ones(nc.mask_rho.shape))
nc['diff_factor'].attrs['long_name'] = "horizontal diffusivity sponge factor" ;
nc['diff_factor'].attrs['valid_min'] = 0. ;
nc['diff_factor'].attrs['coordinates'] = "lon_rho lat_rho" ;

for i in ['visc_factor', 'diff_factor']:
    # for cont in range(30,-2,-1):
    #     nc[i].values[cont-1:-cont+1,-cont] =  4
    #     nc[i].values[cont-1,cont-1:-cont+1] =  4
    #     nc[i].values[-cont,cont-1:-cont+1] =  4

    a = np.ones(nc[i].values.shape)

    cont=1
    for j in range(50,-1,-1):
        a[-j,:] = cont
        a[j,:]  = cont
        a[:,-j] = cont
        a[:, j] = cont

        cont+=0.125
    

    nc[i].values = a

nc['mask_rho'].values[:10,:10] = 0
nc['mask_u'].values[:10,:10] = 0
nc['mask_v'].values[:10,:10] = 0
nc['mask_v'].values[:10,:10] = 0

nc.to_netcdf('/home/dsasaki/Documents/2023/pyroms_tools_newgrid/data/treated_2.012/roms_grid01b_sponge_2.012_smooth1.nc')