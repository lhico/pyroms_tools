from dry_package.plot_schemes import maps
import matplotlib.pyplot as plt
import xarray as xr

nc = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2021/postdoc/202102_roms/pyroms_lhico_tools/data/tpxo_files/grid_tpxo8atlas_30_v1.nc')
nc1 = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2021/postdoc/202102_roms/pyroms_lhico_tools/data/roms_files/sbb_grid_roms.nc')

dxy = 50

fig, ax = plt.subplots(nrows=2,ncols=1)
ax[0].contourf(nc.lon_z[::dxy],  nc.lat_z[::dxy], nc.hz[::dxy, ::dxy].T, cmap=plt.cm.jet)
ax[0].contourf(nc1.lon_rho,  nc1.lat_rho, nc1.h)
ax[1].contourf(nc.nx[::dxy],  nc.ny[::dxy], nc.hz[::dxy, ::dxy].T, cmap=plt.cm.jet)
#
# ax.set_xticklabels(nc.nx.values[::1350])
# ax.set_yticklabels(nc.ny.values[::675])
