import pyroms
from utils import remap
from utils import utils as ut
# from utils import configs
import subprocess
import os, sys
import os.path as osp


os.environ["PYROMS_GRIDID_FILE"] = "gridid.txt"

# -- gets  the information from the config file -- #
# getting the referemce domain from shell 
if len(sys.argv) > 1:
    reference = sys.argv[1]
else:
    reference = 'pbs_202109_glorys'

dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

gridfile      = dicts['grid.grid']
romsgridname  = dicts['gridid.gridname']
output        = dicts['ic.output']
interp_varbs  = dicts['ic.interp_varbs']
interp_varbs2 = dicts['ic.interp_varbs2']
starttime     = dicts['ic.starttime']

# create a gridid.txt (this is a legacy)
# pyroms uses the nc in the variable gridfile
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

# this line uses pyroms gridid.txt
dst_grd = pyroms.grid.get_ROMS_grid(romsgridname)

# -- creates scalar fields -- #
for varb in interp_varbs:
    dfile = f'{dst_grd.name}_{varb}_ic.nc'  # temporary interpolated file
    remap.create_scalar_fields_3D(dst_grd, varb, starttime,
                dst_file=dfile)

# -- creates vector fields -- #
dfileu = f'{dst_grd.name}_u_ic.nc'  # temporary interpolated file name
dfilev = f'{dst_grd.name}_v_ic.nc'  # temporary interpolated file name

# writes intermediate boundary file
remap.create_vector_fields(dst_grd, starttime,dst_fileu=dfileu, dst_filev=dfilev)

# name of the merged initial condition file
ic_file = output
out_file = f'{dst_grd.name}_%s_ic.nc'  # temporary interpolated file


# -- merge files and remove temporary files --#
cont = 0
for istr in interp_varbs2:
    if cont ==0:
        command = ('ncks', '--no_abc', '-O', out_file % istr, ic_file)  # create ic
        cont += 1
    else:
        command = ('ncks', '--no_abc', '-A', out_file % istr, ic_file)  # append
    subprocess.call(command)
    os.remove(out_file % istr)  # remove temp files
print(f'saved to {ic_file}')