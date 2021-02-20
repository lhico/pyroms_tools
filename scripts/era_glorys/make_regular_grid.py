"""
This script creates a simple rectangular roms grid.
"""
import numpy as np
import xarray as xr
import pyroms
from utils import configs
import os.path as osp
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from bathy_smoother import bathy_smoothing, bathy_tools
from utils import utils as ut


# data path
ddir = '../../data/'
odir = '../../data/roms_files'
topofile = '../../data/aux_data/gebco.nc'

rotation_angle = 50  # anticlockwise
# Horizontal grid dimension (these dims are also used in roms config .in file)
Lm = 90
Mm = 40

# vertical coordinate (these coordinates are also used in roms config .in file)
theta_b = 3.0  # change the sigma coordinates close to the bottom
theta_s = 7.0  # change the sigma coordinates close to the surface
Tcline = 250  # parameters
N = 27  # number of depth levels


# # creating a rectangular grid
# lon0=360-50.0 ; lat0=-20.0  # bottom left corners
# lon1=360-50.0 ; lat1=-28.0  # bottom right corner
# lon2=360-35.6 ; lat2=-28.0  # top right corner
# lon3=360-35.6 ; lat3=-20.0  # top left corner

dlon = -2.5
dlat = 0


# creating a rectangular grid
lon0=360-48.0 + dlon ; lat0=-22.0 + dlat  # bottom left corners
lon1=360-48.0 + dlon ; lat1=-29.0 + dlat  # bottom right corner
lon2=360-36.0 + dlon ; lat2=-29.0 + dlat  # top right corner
lon3=360-36.0 + dlon ; lat3=-22.0 + dlat  # top left corner

# define base projection (here mercator)
lon_min = min(lon0, lon1, lon2, lon3)
lon_max = max(lon0, lon1, lon2, lon3)
lon_0 = (lon_min + lon_max) / 2.
lat_min = min(lat0, lat1, lat2, lat3)
lat_max = max(lat0, lat1, lat2, lat3)
lat_0 = (lat_min + lat_max) / 2.

ll = configs.area_params['pcse1']
base = Basemap(projection='merc', llcrnrlon=ll['lonW']+360, llcrnrlat=ll['latS'],
         urcrnrlon=ll['lonE']+360, urcrnrlat=ll['latN'], lat_0=lat_0, lon_0=lon_0,
         resolution='h')

# define the polygon that define the grid boundariees
lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])

 # define vertices rotation directon (one for each vertice 1,0,-1 indicate
 # the direction)
beta = np.array([1, 1, 1, 1])

# rotates the grid
lonrot, latrot = lonp, latp
lonrot, latrot = ut.rotate_coords(lonp, latp, rotation_angle, degrees=True)

#generate the new grid
# Do this if you aren't going to move the grid corners interactively.
hgrd = pyroms.grid.Gridgen(lonrot, latrot, beta, (Mm+3, Lm+3), proj=base)

# use basemap to convert grid coordinate to longitude and latitude
lonv, latv = list(base(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv-1.1, latv-1.3, base)

for xx,yy in base.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy,np.float32)
    vv = np.zeros((xa.shape[0],2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv,mask_value=0)

# Edit the land mask interactively.
#pyroms.grid.edit_mask_mesh(hgrd, proj=base)
#edit_mask_mesh_ij is a faster version using imshow... but no base projection.
coast = pyroms.utility.get_coast_from_map(base)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

# reading a topography file (netcdf)
nc = xr.open_dataset(topofile)
topo = nc['elevation'][::5,::5].values
lons = nc['lon'][::5].values
lats = nc['lat'][::5].values

lonm, latm = np.meshgrid(lons, lats)
lonm = lonm.ravel()+360
latm = latm.ravel()

# depth positive
topo = -topo

# fix minimum depth
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
h = griddata((lonm.ravel(),latm.ravel()),topo.ravel(),
             (hgrd.lon_rho,hgrd.lat_rho), method='linear')

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# save raw bathymetry
hraw = h.copy()

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 1.0
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

vgrd = pyroms.vgrid.s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw)
# print((vgrd.Cs_r*4000).astype('int'))
# plt.plot(vgrd.s_rho, vgrd.Cs_r, marker='x')

# ROMS grid
grd_name = 'SBB4'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
outfile = osp.join(odir, 'sbb_grid_roms.nc')
pyroms.grid.write_ROMS_grid(grd, filename=outfile)
print(f'file saved at {outfile}')
