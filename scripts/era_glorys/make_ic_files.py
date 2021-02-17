import pyroms
from utils import CGrid_glorys
from utils import remap
from utils import configs
import subprocess
import os
import os.path as osp
from utils import configs

datapath = '../../data'

name = 'GLORYS_SWATL'
area='regional'

# -- the xrange and yrange MUST be the same of make_remap_weights_file.py -- #
xrange=configs.area_params['pcse1']['xrange']
yrange=configs.area_params['pcse1']['yrange']

# -- files used to write the initial conditions -- #
inifile = 'phy-001-024-daily_55W30W35S15S_2019-08-01_12:00:00.nc'  # data
mskfile = 'glo_mask_bathy2.nc'  # mask dataset (created by make_ref_mask.py)

# -- the xrange and yrange MUST be the same of make_remap_weights_file.py -- #
xrange = configs.area_params['pcse1']['xrange']
yrange = configs.area_params['pcse1']['yrange']

file = osp.join(datapath, 'ocean_model', inifile)
fpath = osp.join(datapath, 'aux_data', mskfile)
outdir = osp.join(datapath, 'roms_files')

starttime = '2019-08-01T12:00:00'

wfiles = [
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_t_to_rho.nc',
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_u_to_rho.nc',
    'remap_weights_GLORYS_SWATL_to_SBB4_bilinear_v_to_rho.nc'
    ]

# -- cut glorys 030 mask file to the same coords of the downloaded dataset -- #
# extent = configs.area_params['pcse1']
# extent = configs.area_params['wSAtl0']
# # fpath = '/home/otel/Desktop/GLO-MFC_001_030_mask_bathy.nc'
# nc = xr.open_dataset(fpath)
# ncut = nc.sel(longitude=slice(extent['lonW'], extent['lonE']),
#               latitude=slice(extent['latS'], extent['latN']))
# ncut.to_netcdf('glo_mask_bathy.nc')

# -- gets the source and destiny grids --#
src_grd = CGrid_glorys.A2CGrid(fpath, xrange=xrange, yrange=yrange, lonoffset=360)
dst_grd = pyroms.grid.get_ROMS_grid('SBB4')

# -- interpolates scalar fields -- #
for varb in ['zos', 'so', 'thetao']:
    dfile = f'{dst_grd.name}_{varb}_ic.nc'  # temporary interpolated file
    remap.remap(file, src_grd, dst_grd, varb, starttime, dst_file=dfile,
                weight_file=wfiles[0])

# -- interpolates vector fields -- #
dfileu = f'{dst_grd.name}_u_ic.nc'  # temporary interpolated file name
dfilev = f'{dst_grd.name}_v_ic.nc'  # temporary interpolated file name

# interpolating
remap.remap_uv(file, src_grd, dst_grd, starttime,
               dst_fileu=dfileu, dst_filev=dfilev, wts_file=wfiles)

# name of the merged initial condition file
ic_file = bdry_file = osp.join(outdir, f'{dst_grd.name}_ic.nc')
out_file = f'{dst_grd.name}_%s_ic.nc'  # temporary interpolated file

# -- merge files and remove temporary files --#
cont = 0
for istr in ['zos', 'so', 'thetao', 'u', 'v']:
    if cont ==0:
        command = ('ncks', '--no_alphabetize', '-O', out_file % istr, ic_file)  # create ic
        cont += 1
    else:
        command = ('ncks', '--no_alphabetize', '-A', out_file % istr, ic_file)  # append
    subprocess.call(command)
    os.remove(out_file % istr)  # remove temp files
