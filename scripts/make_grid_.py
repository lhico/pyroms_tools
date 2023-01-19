import sys
import matplotlib.pyplot as plt
import os.path as osp
import numpy as np
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import xarray as xr
# from utils import configs
import pyroms
from scipy import signal
from utils import utils as ut

def hgrid(lon_rho, lat_rho):
    # Grid positions
    lonv = lon_rho.copy()
    latv = lat_rho.copy()

    print(lonv.shape)
    # basemap 
    base = Basemap(projection='merc',
                llcrnrlon=lonv.min()-2, llcrnrlat=latv.min()-2,
                urcrnrlon=lonv.max()+2, urcrnrlat=latv.max()+2, lat_0=latv.mean(), lon_0=lonv.mean(),
            resolution='h')

    # horizontal grid object
    hgrd = pyroms.grid.CGrid_geo(lonv, latv, base)

    # generate the mask based on the coastlines
    for xx,yy in base.coastpolygons:
        xa = np.array(xx, np.float32)
        ya = np.array(yy,np.float32)
        vv = np.zeros((xa.shape[0],2))
        vv[:, 0] = xa
        vv[:, 1] = ya
        hgrd.mask_polygon(vv,mask_value=0)
    return hgrd, base

def h_bathymetry(topo, lon, lat, hgrd):
    # -- prepare bathymetry -- #
    hmin = 5
    # topo = nc.hraw.values
    topo = np.where(topo<hmin, hmin, topo)
    # lon = nc.lon_rho.values
    # lat = nc.lat_rho.values

    # interpolate new bathymetry
    h = griddata((lon.ravel(),lat.ravel()),topo.ravel(),(hgrd.lon_rho,hgrd.lat_rho), method='linear')


    # ensure that depth is always deeper than hmin
    h = np.where(h < hmin, hmin, h)

    idx = np.where(hgrd.mask_rho == 0)
    h[idx] = hmin

    # save raw bathymetry
    hraw = h.copy()

    # -- check bathymetry -- #
    # check bathymetry roughness
    RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
    print('Max Roughness value is: ', RoughMat.max())

    # # # smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
    # rx0_max = 0.4
    # h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

    # check bathymetry roughness again
    RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
    print('Max Roughness value is: ', RoughMat.max())

    return h,hraw

def rotate_coords(xm, ym, ang_rot, degrees=False):
    xm_mean = np.mean(xm)
    ym_mean = np.mean(ym)

    xm1 = xm.copy() - xm_mean
    ym1 = ym.copy() - ym_mean

    ang_rot_rad = ang_rot
    if degrees:
        ang_rot_rad = ang_rot/180 *np.pi

    xrot = xm1*np.cos(ang_rot_rad) - ym1*np.sin(ang_rot_rad)
    yrot = xm1*np.sin(ang_rot_rad) + ym1*np.cos(ang_rot_rad)

    xrot += xm_mean
    yrot += ym_mean
    return xrot, yrot

def interpolate_bathymetry(ddir, hgrd):
    window = np.ones([5,5])
    window = window/window.sum()
    nc = xr.open_dataset(ddir)
    topo = nc['elevation'].values #  np.loadtxt(os.path.join(ddir, 'etopo20data.gz'))
    topo1 = signal.convolve2d(topo, window, 'same')
    topo = topo1[::5,::5]
    lons = nc['lon'][::5].values  # np.loadtxt(os.path.join(ddir, 'etopo20lons.gz'))
    lats = nc['lat'][::5].values  # np.loadtxt(os.path.join(ddir, 'etopo20lats.gz'))

    lonm, latm = np.meshgrid(lons,lats)
    lonm = lonm.ravel()
    latm = latm.ravel()

    # depth positive
    topo = -topo

    # fix minimum depth
    hmin = 5
    topo = np.where(topo < hmin, hmin, topo)

    # interpolate new bathymetry
    lon, lat = np.meshgrid(lons, lats)
    h = griddata((lonm.ravel(),latm.ravel()),topo.ravel(),(hgrd.lon,hgrd.lat), method='linear')
    return h

if __name__ == '__main__':
    # plt.close('all')

    # -- gets  the information from the config file -- #
    # getting the referemce domain from shell 
    if len(sys.argv) > 1:
        reference = sys.argv[1]
    else:
        reference = 'pbs_202109_glorys'

    dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
    dicts = dicts[reference]
    
    # -- paths -- #s
    odir  = dicts['output_dir'] 
    bfile = dicts['bathy_file']

    # -- horizontal grid parameters -- #
    dxdy        = dicts['grid.dxdy']
    x0,x1,y0,y1 = dicts['grid.WESN']
    xoffset     = dicts['grid.xoffset']
    yoffset     = dicts['grid.yoffset']
    rot         = dicts['grid.rot']

    # -- vertical grid parameters -- #
    N       = dicts['grid.N'] 
    theta_s = dicts['grid.theta_s']
    theta_b = dicts['grid.theta_b']
    Tcline  = dicts['grid.Tcline']


    #  -- grid rotation -- #
    # x = np.arange(-44,-40,1)+1
    x = np.arange(x0,x1,dxdy) + xoffset
    y = np.arange(y0,y1,dxdy) + yoffset


    xm, ym = np.meshgrid(x,y)
    xrot, yrot =rotate_coords(xm, ym, rot)
    hgrd, base = hgrid(xrot, yrot)

    #pyroms.grid.edit_mask_mesh(hgrd, proj=base)
    #edit_mask_mesh_ij is a faster version using imshow... but no base projection.
    coast = pyroms.utility.get_coast_from_map(base)
    pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

    topo = interpolate_bathymetry(bfile, hgrd)

    # topo = -hraw1    # getting bathymetry (in this case the bathymetry was in the hromsg)
    lon = xrot  # getting position
    lat = yrot  # getting position
    h,hraw = h_bathymetry(topo, lon, lat, hgrd)  # bathymetry [object (?)]

    plt.close('all')
    # vertical coordinate
    vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
    # print((vgrd.Cs_r*4000).astype('int'))
    # plt.plot(vgrd.s_rho, vgrd.Cs_r, marker='x')


    # ROMS grid
    grd_name = 'BaciaSantos'
    grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

    # write grid to netcdf file
    pyroms.grid.write_ROMS_grid(grd, filename=osp.join(odir, 'roms_grid00.nc'))


