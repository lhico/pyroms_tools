# this script intends to merge more than one file for 
# each variable from ERA5. This happens particularly
# when dealing with a period to simulate that covers
# more than one year
import xarray as xr
from glob import glob
import numpy as np

datadir = '/home/danilo/data/reanalysis/era5'
experim = 'deproas_v2002'
outpPat = '_'
numberF = 2

nfiles = sorted(glob(f"{datadir}/*.nc"))

# taking files every numberF
old = 0
for i in np.arange(numberF, len(nfiles), numberF):
    _tmp = nfiles[old:i]
    old = i

    ds = xr.open_mfdataset(_tmp, decode_cf=False)

    # new name
    outp = _tmp[0][:_tmp[0].rfind("_")] + '.nc'
    # saving new file
    ds.to_netcdf(outp)
