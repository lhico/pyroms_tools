import pyroms
from utils import CGrid_glorys
from utils import remap_bdry
from utils import configs
import subprocess
import os
import os.path as osp
import pandas as pd
import glob

datapath = '../../data'

name = 'GLORYS_SWATL'
area='regional'
startTime = '2019-08-01T12:00:00'
endTime   = '2019-09-01T12:00:00'
freq='1D'

# -- files used to write the initial conditions -- #
inifile = '*.nc'  # data
mskfile = 'glo_mask_bathy2.nc'  # mask dataset (created by make_ref_mask.py)

# -- the xrange and yrange MUST be the same of make_remap_weights_file.py -- #
xrange = configs.area_params['pcse1']['xrange']
yrange = configs.area_params['pcse1']['yrange']

file = sorted(glob.glob(osp.join(datapath, 'ocean_model', inifile)))
fpath = osp.join(datapath, 'aux_data', mskfile)
outdir = osp.join(datapath, 'roms_files')

timerange = pd.date_range(startTime, endTime, freq=freq)


wfiles = [
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_t_to_rho.nc',
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_u_to_rho.nc',
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_v_to_rho.nc'
    ]

src_grd = CGrid_glorys.A2CGrid(fpath, xrange=xrange, yrange=yrange, lonoffset=360)
dst_grd = pyroms.grid.get_ROMS_grid('SBB4')  # read data from the gridid.txt


for it, tfile in zip(timerange, file):
    itstr = str(it).replace(' ', 'T')

    # -- interpolates scalar fields -- #
    for varb in ['zos', 'so', 'thetao']:
        dfile = f'{dst_grd.name}_{varb}_bdry.nc'  # temporary interpolated file
        remap_bdry.remap_bdry(tfile, src_grd, dst_grd, varb, itstr,
                              dst_file=dfile, weight_file=wfiles[0])

    # -- interpolates vector fields -- #
    dfileu = f'{dst_grd.name}_u_bdry.nc'  # temporary interpolated file name
    dfilev = f'{dst_grd.name}_v_bdry.nc'  # temporary interpolated file name

    # interpolating
    remap_bdry.remap_bdry_uv(tfile, src_grd, dst_grd, itstr,
                   dst_fileu=dfileu, dst_filev=dfilev, wts_file=wfiles)

    # name of the merged initial condition file
    bdry_file = osp.join(outdir, f'{dst_grd.name}_bdry_{itstr}.nc')
    out_file = f'{dst_grd.name}_%s_bdry.nc'  # temporary interpolated file

    # -- merge files and remove temporary files --#
    cont = 0
    for istr in ['zos', 'so', 'thetao', 'u', 'v']:
        if cont ==0:
            command = ('ncks', '--no_alphabetize', '-O', out_file % istr, bdry_file)  # create ic
            cont += 1
        else:
            command = ('ncks', '--no_alphabetize', '-A', out_file % istr, bdry_file)  # append
        subprocess.call(command)
        os.remove(out_file % istr)  # remove temp files
