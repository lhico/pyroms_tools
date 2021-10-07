import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import xesmf as xe

odir = '/home/otel/Dropbox/trabalho_irado/2021/postdoc/202101_caracterizacao_ambiental_PCSE/roms_tools_projetoBS1.2/data/roms_files/ecom/20210915_clim/output/ecom_glorys_ic.nc'

ncecom = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2021/postdoc/2021_data/ecom/gcmplt.nc')
ecomgrid = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2021/postdoc/2021_data/ecom/grid/ecom_modelgrid.nc')
ncecom.lon.values[1:-1,1:-1] = ecomgrid.lon.values
ncecom.lat.values[1:-1,1:-1] = ecomgrid.lat.values
ncecom = ncecom.drop(['xpos', 'ypos'])
ncecom = ncecom.rename({'xpos':'x', 'ypos':'y'})
ncecom = ncecom.isel(x=slice(1,-1), y=slice(1,-1))

ncroms = xr.open_dataset(odir)
ncroms0 = ncroms.copy()
ncroms = ncroms.rename({'lon_rho':'lon',
                        'lat_rho':'lat',
                        'xi_rho': 'x',
                        'eta_rho': 'y',
                        'ocean_time': 'time',
                        's_rho': 'sigma'})  # rename variables so xesmf understand them

for varb in ['salt', 'temp']:
    regridder = xe.Regridder(ncecom[varb], ncroms, 'bilinear', extrap_method='nearest_s2d')
    interpolated = regridder(ncecom[varb][0])  # interpolating

    aux = np.flipud(interpolated.values.copy())
    aux1 = ncroms0[varb].values[0].copy()

    aux1[~np.isnan(aux)] = aux[~np.isnan(aux)]


    ncroms0[varb].values[0] = aux1

for mask in ['rho', 'u', 'v']:
    ncroms0[f'mask_{mask}'].values[0,-1] = 0
    ncroms0[f'mask_{mask}'].values[-1,-1] = 0

ncroms0['u'].values[:] = 0
ncroms0['v'].values[:] = 0
ncroms0['ubar'].values[:] = 0
ncroms0['vbar'].values[:] = 0


ncroms0.to_netcdf('ecom_lhico.nc')