import pyroms
from utils import CGrid_glorys
from utils import remap_bdry
import subprocess
import os
import os.path as osp
import pandas as pd
import glob
import sys
from utils import utils as ut
os.environ["PYROMS_GRIDID_FILE"] = "gridid.txt"

# -- gets  the information from the config file -- #
# getting the referemce domain from shell 
if len(sys.argv) > 1:
    reference = sys.argv[1]
else:
    reference = 'pbs_202109_glorys'

dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

romsgridname  = dicts['gridid.gridname']       # roms grid name
startTime     = dicts['bndry.dummytime']       # dummy time 
outdir        = dicts['output_dir']            # directory where we save stuff
gridfile      = dicts['grid.grid']

# create a gridid.txt (this is a legacy from the original pyroms)
ut.create_grid_file(
    romsgridname,
    romsgridname,
    gridfile,
    dicts['grid.N'],
    'roms',
    2,
    dicts['grid.theta_s'],
    dicts['grid.theta_b'],
    dicts['grid.Tcline']
)


# -- files used to write the initial conditions -- #
timerange = pd.date_range(startTime, startTime, freq='1D')

dst_grd = pyroms.grid.get_ROMS_grid(romsgridname)  # read data from the gridid.txt

# using only one file to create an empty boundary grid template, which will be filled with
# global outputs by the scripts_xesmf/make_bdry_file_replace.py
it = timerange[0]
itstr = str(it).replace(' ', 'T')

# -- interpolates vector fields -- #
dfileu = f'{dst_grd.name}_u_bdry.nc'  # temporary interpolated file name
dfilev = f'{dst_grd.name}_v_bdry.nc'  # temporary interpolated file name

# interpolating
remap_bdry.make_bdry_uv_file(dst_grd, itstr,dst_fileu=dfileu, dst_filev=dfilev)

# -- interpolates scalar fields -- #
for varb in ['so', 'thetao','zos']:
    dfile = f'{dst_grd.name}_{varb}_bdry.nc'  # temporary interpolated file
    remap_bdry.make_bdry_rho_file(dst_grd, varb, itstr,dst_file=dfile)

# name of the merged initial condition file
bdry_file = osp.join(outdir, 'empty_bdry.nc')
out_file = f'{dst_grd.name}_%s_bdry.nc'  # temporary interpolated file

# -- merge files and remove temporary files --#
cont = 0
for istr in ['zos', 'so', 'thetao', 'u', 'v']:
    if cont ==0:
        command = ('ncks', '--no_abc', '-O', out_file % istr, bdry_file)  # create ic
        cont += 1
    else:
        command = ('ncks', '--no_abc', '-A', out_file % istr, bdry_file)  # append
    subprocess.call(command)
    os.remove(out_file % istr)  # remove temp files

