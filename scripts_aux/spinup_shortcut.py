# roms tends to produce  a lot of inertial waves before reaching
# the equilibrium. A way I found to shortcut this is to let the
# model to run from resting state and after a few inertial periods
# take the median in time of u and v (this is not u_eastward, nor
# v_eastward). In the experinces with western boundary currents
# this result made the results converge much quicker.
import xarray as xr
import numpy as np

# ffile = 'roms_avg.nc'
# icfile = 'input2/nestcoarse_glorys_ic2.nc'
# outfile = 'input2/nestcoarse_glorys_ic3.nc'

# ffile = 'roms_avg_nest.nc'
# icfile = 'input2/nestref1_glorys_ic2.nc'
# outfile = 'input2/nestref1_glorys_ic3.nc'

# ffile = 'roms_avg_nest2.nc'
# icfile = 'input2/nestref2_glorys_ic2.nc'
# outfile = 'input2/nestref2_glorys_ic3.nc'

# ffile = 'roms_avg.nc'
# icfile = 'input/pBS_glorys_ic.nc'
# outfile = 'input/pBS_glorys_ic_spup2.nc'

ffile = 'roms_avg.nc'
icfile = 'input/pbs_202109_smooth_ic2.nc'
outfile = 'input/pbs_202109_smooth_ic3_spup.nc'


# ffile = 'roms_avg.nc'
# icfile = 'input/woa_ic.nc'
#  outfile = 'input/woa_ic_spup.nc'


ncavg = xr.open_dataset(ffile)
nc_ic = xr.open_dataset(icfile)

u = ncavg.u[:20]
v = ncavg.v[:20]

nc_ic.u.values = [np.nanmedian(u, axis=0)]
nc_ic.v.values = [np.nanmedian(v, axis=0)]

nc_ic.to_netcdf(outfile)