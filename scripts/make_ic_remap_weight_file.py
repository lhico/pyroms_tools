import pyroms
import pyroms_toolbox
from utils import CGrid_glorys
import os.path as osp
import xarray as xr
import numpy as np
from utils import utils as ut
import os
os.environ["PYROMS_GRIDID_FILE"] = "/home/lhico/pyroms_tools/configs/gridid.txt"



# -- gets  the information from the config file -- #
reference = 'pbs_202109_glorys'
dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

mpath  = dicts['ic.grid_file']
romsgridname   = dicts['ic.gridname']
srcgridname = dicts['ic.srcname']
xrange = dicts['ic.source_xrange']
yrange = dicts['ic.source_yrange']
varstr = dicts['ic.varstr']

srcgrd = CGrid_glorys.A2CGrid(mpath, xrange=xrange, yrange=yrange, name=srcgridname, **varstr)

dstgrd = pyroms.grid.get_ROMS_grid(romsgridname)

pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='t')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='u')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{srcgridname}_t.nc'
grid2_file = f'remap_grid_{romsgridname}_rho.nc'
interp_file1 = f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_rho.nc'
interp_file2 = f'remap_weights_{romsgridname}_to_{srcgridname}_bilinear_rho_to_t.nc'
map1_name = f'{srcgridname} to {romsgridname} Bilinear Mapping'
map2_name = f'{romsgridname} to {srcgridname} Bilinear Mapping'
num_maps = 1
map_method = 'distwgt'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{srcgridname}_u.nc'
grid2_file = f'remap_grid_{romsgridname}_rho.nc'
interp_file1 = f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_u_to_rho.nc'
interp_file2 = f'remap_weights_{romsgridname}_to_{srcgridname}_bilinear_rho_to_u.nc'
map1_name = f'{srcgridname} to {romsgridname} Bilinear Mapping'
map2_name = f'{romsgridname} to {srcgridname} Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{srcgridname}_v.nc'
grid2_file = f'remap_grid_{romsgridname}_rho.nc'
interp_file1 = f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_v_to_rho.nc'
interp_file2 = f'remap_weights_{romsgridname}_to_{srcgridname}_bilinear_rho_to_v.nc'
map1_name = f'{srcgridname} to {romsgridname} Bilinear Mapping'
map2_name = f'{romsgridname} to {srcgridname} Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{srcgridname}_t.nc'
grid2_file = f'remap_grid_{romsgridname}_u.nc'
interp_file1 = f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_u.nc'
interp_file2 = f'remap_weights_{romsgridname}_to_{srcgridname}YS_bilinear_u_to_t.nc'
map1_name = f'{srcgridname} to {romsgridname} Bilinear Mapping'
map2_name = f'{romsgridname} to {srcgridname} Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{srcgridname}_t.nc'
grid2_file = f'remap_grid_{romsgridname}_v.nc'
interp_file1 = f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_v.nc'
interp_file2 = f'remap_weights_{romsgridname}_to_{srcgridname}YS_bilinear_v_to_t.nc'
map1_name = f'{srcgridname} to {romsgridname} Bilinear Mapping'
map2_name = f'{romsgridname} to {srcgridname} Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')
