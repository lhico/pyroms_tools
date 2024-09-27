from utils import tpxo8 as tpxtoroms
from utils import utils as ut
import sys
import datetime

# -- gets  the information from the config file -- #
# getting the referemce domain from shell 
if len(sys.argv) > 1:
    reference = sys.argv[1]
else:
    reference = 'pbs_202109_glorys'

dicts = ut._get_dict_paths('../configs/grid_config_esmf.txt')
dicts = dicts[reference]

harmonic_names = dicts['tpxo.harmonic_names']
tpxofile_data  = dicts['tpxo.dict']
t0=datetime.datetime(*dicts['tpxo.datetime'])
ndays = dicts['tpxo.ndays']
# 
ncROMSgrdname = dicts['grid_dir']
ncoutFilename = dicts['tpxo.outfile']
tpxtoroms.tpxo2roms(t0, ncROMSgrdname,
                   harmonic_names,
                   tpxofile_data,
                   ncoutFilename=ncoutFilename,
                   ndays=ndays)