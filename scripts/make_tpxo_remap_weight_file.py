import pyroms
import pyroms_toolbox
import CGrid_TPXO8
from utils import utils as ut
import os
os.environ["PYROMS_GRIDID_FILE"] = "/home/lhico/pyroms_tools/configs/gridid.txt"


# -- gets  the information from the config file -- #
reference = 'pbs_202109_glorys'
dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

gridname = dicts['ic.gridname']
tpxogrid = dicts['tpxo.grid']
xrange   = dicts['tpxo.xrange']
yrange   = dicts['tpxo.yrange']

srcgrd = CGrid_TPXO8.get_nc_CGrid_TPXO8(tpxogrid, xrange=xrange, yrange=yrange, lonoffset=-360)
#srcgrd_lr = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6', \
#      xrange=(920, 1020), yrange=(1030, 1300))
dstgrd = pyroms.grid.get_ROMS_grid(gridname)

# make remap grid file for scrip
CGrid_TPXO8.make_remap_grid_file(srcgrd, Cpos='t')
CGrid_TPXO8.make_remap_grid_file(srcgrd, Cpos='u')
CGrid_TPXO8.make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_t.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_t.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


grid1_file = 'remap_grid_' + srcgrd.name + '_u.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_u_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_u.nc'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)


grid1_file = 'remap_grid_' + srcgrd.name + '_v.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_v_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_v.nc'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method)
