import matplotlib
# matplotlib.use('Agg')
import pyroms
import pyroms_toolbox
from utils import CGrid_glorys
import os.path as osp
from utils import configs
#import CGrid_GLORYS

datapath = '../../data'
fpath = osp.join(datapath, 'aux_data','glo_mask_bathy2.nc')
xrange=configs.area_params['pcse1']['xrange']
yrange=configs.area_params['pcse1']['yrange']

name = 'GLORYS_SWATL'
outgrid = 'SBB4'
# load the grid
# srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_GLORYS('/archive/u1/uaf/kate/GLORYS/GL2V1_mesh_mask_new.nc', name='GLORYS_ARCTIC2', area='npolar', ystart=690)
srcgrd = CGrid_glorys.A2CGrid(fpath, xrange=xrange, yrange=yrange, name=name, lonoffset=360)
dstgrd = pyroms.grid.get_ROMS_grid('SBB4')

# make remap grid file for scrip
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='t')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='u')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{name}_t.nc'
grid2_file = f'remap_grid_{outgrid}_rho.nc'
interp_file1 = f'remap_weights_{name}_to_{outgrid}_bilinear_t_to_rho.nc'
interp_file2 = f'remap_weights_{outgrid}_to_{name}_bilinear_rho_to_t.nc'
map1_name = f'GLORYS to {outgrid} Bilinear Mapping'
map2_name = f'{outgrid} to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{name}_u.nc'
grid2_file = f'remap_grid_{outgrid}_rho.nc'
interp_file1 = f'remap_weights_{name}_to_{outgrid}_bilinear_u_to_rho.nc'
interp_file2 = f'remap_weights_{outgrid}_to_{name}_bilinear_rho_to_u.nc'
map1_name = f'GLORYS to {outgrid} Bilinear Mapping'
map2_name = f'{outgrid} to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{name}_v.nc'
grid2_file = f'remap_grid_{outgrid}_rho.nc'
interp_file1 = f'remap_weights_{name}_to_{outgrid}_bilinear_v_to_rho.nc'
interp_file2 = f'remap_weights_{outgrid}_to_{name}_bilinear_rho_to_v.nc'
map1_name = f'GLORYS to {outgrid} Bilinear Mapping'
map2_name = f'{outgrid} to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{name}_t.nc'
grid2_file = f'remap_grid_{outgrid}_u.nc'
interp_file1 = f'remap_weights_{name}_to_{outgrid}_bilinear_t_to_u.nc'
interp_file2 = f'remap_weights_{outgrid}_to_{name}YS_bilinear_u_to_t.nc'
map1_name = f'GLORYS to {outgrid} Bilinear Mapping'
map2_name = f'{outgrid} to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')


# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = f'remap_grid_{name}_t.nc'
grid2_file = f'remap_grid_{outgrid}_v.nc'
interp_file1 = f'remap_weights_{name}_to_{outgrid}_bilinear_t_to_v.nc'
interp_file2 = f'remap_weights_{outgrid}_to_{name}YS_bilinear_v_to_t.nc'
map1_name = f'GLORYS to {outgrid} Bilinear Mapping'
map2_name = f'{outgrid} to GLORYS Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file, \
              interp_file1, interp_file2, map1_name, \
              map2_name, num_maps, map_method, \
              grid1_periodic='.false.', grid2_periodic='.false.')
