from utils import utils as ut
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
output        = dicts['nudg.output']

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

onudg = ut.nudgcoef(romsgridname)


# strong restoring
east   = dicts['nudg.east']
west   = dicts['nudg.west']
north  = dicts['nudg.north']
south  = dicts['nudg.south']

tracer_timescales = dicts['nudg.tracertscale']
onudg(east,west,north,south,tracer_timescales,foutname=output)
