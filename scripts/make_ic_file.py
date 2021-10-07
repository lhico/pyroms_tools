import pyroms
from utils import CGrid_glorys
from utils import remap
from utils import utils as ut
# from utils import configs
import subprocess
import os
import os.path as osp
import os
os.environ["PYROMS_GRIDID_FILE"] = "/home/otel/Dropbox/trabalho_irado/2021/postdoc/202101_caracterizacao_ambiental_PCSE/roms_tools_projetoBS1.2/configs/gridid.txt"

# -- gets  the information from the config file -- #
reference = 'pbs_202109_glorys'
dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
dicts = dicts[reference]

icgridfile   = dicts['ic.grid_file'] 
icfile       = dicts['ic.input_file']
outdir       = dicts['output_dir']
romsgridname = dicts['ic.gridname']
srcgridname  = dicts['ic.srcname']
interp_varbs  = dicts['ic.interp_varbs']
interp_varbs2 = dicts['ic.interp_varbs2']
varstr       = dicts['ic.varstr']
xrange       = dicts['ic.source_xrange']
yrange       = dicts['ic.source_yrange']
scalar_dict  = dicts['scalar_dict']
starttime    = dicts['ic.starttime']

# area='regional'

# -- the xrange and yrange MUST be the same of make_remap_weights_file.py -- #

# -- files used to write the initial conditions -- #



wfiles = [
    f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_rho.nc',
    f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_rho.nc',
    f'remap_weights_{srcgridname}_to_{romsgridname}_bilinear_t_to_rho.nc'
    ]


# -- gets the source and destiny grids --#    nc.to_netcdf('woa_clim_extend.nc')
src_grd = CGrid_glorys.A2CGrid(icgridfile, xrange=xrange, yrange=yrange, lonoffset=0, **varstr)
dst_grd = pyroms.grid.get_ROMS_grid(romsgridname)

# -- interpolates scalar fields -- #
for varb in interp_varbs:
    dfile = f'{dst_grd.name}_{varb}_ic.nc'  # temporary interpolated file
    remap.remap(icfile, src_grd, dst_grd, varb, starttime, scalar_dict,
                dst_file=dfile,
                weight_file=wfiles[0])

# -- interpolates vector fields -- #
dfileu = f'{dst_grd.name}_u_ic.nc'  # temporary interpolated file name
dfilev = f'{dst_grd.name}_v_ic.nc'  # temporary interpolated file name

# interpolating
remap.remap_uv(icfile, src_grd, dst_grd, starttime,
               dst_fileu=dfileu, dst_filev=dfilev, wts_file=wfiles)

# name of the merged initial condition file
ic_file = bdry_file = osp.join(outdir, f'{dst_grd.name}_{srcgridname}_ic.nc')
out_file = f'{dst_grd.name}_%s_ic.nc'  # temporary interpolated file


# -- merge files and remove temporary files --#
cont = 0
for istr in interp_varbs2:
    if cont ==0:
        command = ('ncks', '-a', '-O', out_file % istr, ic_file)  # create ic
        cont += 1
    else:
        command = ('ncks', '-a', '-A', out_file % istr, ic_file)  # append
    subprocess.call(command)
    os.remove(out_file % istr)  # remove temp files
print(f'saved to {ic_file}')