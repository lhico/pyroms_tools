import pyroms
from utils import CGrid_glorys
from utils import remap_bdry
import subprocess
import os
import os.path as osp
import pandas as pd
import glob
from utils import utils as ut
os.environ["PYROMS_GRIDID_FILE"] = "/home/lhico/pyroms_tools/configs/gridid.txt"

# plt.close('all')
reference = 'swatl_2022'
dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]


datapath  = dicts['bndry.srcfiles']       # dataset source  (string)
srcname   = dicts['bndry.srcname']        # name of the dataset source (string)
gridname  = dicts['bndry.gridname']       # roms grid name
startTime = dicts['bndry.daterange'][0]   # string for pd.date_range
endTime   = dicts['bndry.daterange'][1]   # string for pd.date_range
freq      = dicts['bndry.daterange'][2]   # string for pd.date_range
varstr    = dicts['bndry.varstr']         # coordinates names of the source dataset
outdir    = dicts['output_dir']           # directory where we save stuff

# -- files used to write the initial conditions -- #
mskfile = dicts['bndry.srcgrid']  # mask dataset (created by make_ref_mask.py)

# -- the xrange and yrange MUST be the same of make_remap_weights_file.py -- #
xrange = dicts['bndry.xrange']  # x indices limits of the source dataset
yrange = dicts['bndry.yrange']  # x indices limits of the source dataset

file = sorted(glob.glob(datapath))
fpath = osp.join(datapath, mskfile)

timerange = pd.date_range(startTime, endTime, freq=freq)


wfiles = [
    f'remap_weights_{srcname}_to_{gridname}_bilinear_t_to_rho.nc',
    f'remap_weights_{srcname}_to_{gridname}_bilinear_u_to_rho.nc',
    f'remap_weights_{srcname}_to_{gridname}_bilinear_v_to_rho.nc'
    ]

src_grd = CGrid_glorys.A2CGrid(fpath, xrange=xrange, yrange=yrange, lonoffset=0, **varstr)
dst_grd = pyroms.grid.get_ROMS_grid(gridname)  # read data from the gridid.txt


for it, tfile in zip(timerange, file):
    itstr = str(it).replace(' ', 'T')

    # -- interpolates vector fields -- #
    dfileu = f'{dst_grd.name}_u_bdry.nc'  # temporary interpolated file name
    dfilev = f'{dst_grd.name}_v_bdry.nc'  # temporary interpolated file name

    # interpolating
    remap_bdry.make_bdry_uv_file(tfile, src_grd, dst_grd, itstr,
                   dst_fileu=dfileu, dst_filev=dfilev, wts_file=wfiles)


    # -- interpolates scalar fields -- #
    for varb in ['so', 'thetao','zos']:
        dfile = f'{dst_grd.name}_{varb}_bdry.nc'  # temporary interpolated file
        remap_bdry.make_bdry_rho_file(tfile, src_grd, dst_grd, varb, itstr,
                              dst_file=dfile, weight_file=wfiles[0])


    # name of the merged initial condition file
    bdry_file = osp.join(outdir, f'{dst_grd.name}_bdry_{itstr}.nc')
    out_file = f'{dst_grd.name}_%s_bdry.nc'  # temporary interpolated file

    # -- merge files and remove temporary files --#
    cont = 0
    for istr in ['zos', 'so', 'thetao', 'u', 'v']:
        if cont ==0:
            command = ('ncks', '-a', '-O', out_file % istr, bdry_file)  # create ic
            cont += 1
        else:
            command = ('ncks', '-a', '-A', out_file % istr, bdry_file)  # append
        subprocess.call(command)
        os.remove(out_file % istr)  # remove temp files
