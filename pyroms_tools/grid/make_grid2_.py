# import matplotlib.pyplot as plt
# import os.path as osp
# # from pyroms import _iso
# import numpy as np
# from mpl_toolkits.basemap import Basemap, shiftgrid
# from scipy.interpolate import griddata
# import xarray as xr
# # from utils import configs
# import pyroms
# from scipy import ndimage
# from scipy import signal
# from bathy_smoother import bathy_tools, bathy_smoothing
# from utils import utils as ut
# import sys
# import pyproj


# def hgrid(lon_rho, lat_rho):
#     # Grid positions
#     lonv = lon_rho.copy()
#     latv = lat_rho.copy()

#     print(lonv.shape)
#     # basemap 
#     base = Basemap(projection='aeqd',
#                 llcrnrlon=lonv.min()-2, llcrnrlat=latv.min()-2,
#                 urcrnrlon=lonv.max()+2, urcrnrlat=latv.max()+2, lat_0=latv.mean(), lon_0=lonv.mean(),
#             resolution='h')

#     # horizontal grid object
#     hgrd = pyroms.grid.CGrid_geo(lonv, latv, base)

#     # generate the mask based on the coastlines
#     for xx,yy in base.coastpolygons:
#         xa = np.array(xx, np.float32)
#         ya = np.array(yy,np.float32)
#         vv = np.zeros((xa.shape[0],2))
#         vv[:, 0] = xa
#         vv[:, 1] = ya
#         hgrd.mask_polygon(vv,mask_value=0)
#     return hgrd, base


# def h_bathymetry(topo, lon, lat, hgrd):
#     # -- prepare bathymetry -- #
#     hmin = 5
#     # topo = nc.hraw.values
#     topo = np.where(topo<hmin, hmin, topo)
#     # lon = nc.lon_rho.values
#     # lat = nc.lat_rho.values

#     # interpolate new bathymetry
#     h = griddata((lon.ravel(),lat.ravel()),topo.ravel(),(hgrd.lon_rho,hgrd.lat_rho), method='linear')


#     # ensure that depth is always deeper than hmin
#     h = np.where(h < hmin, hmin, h)

#     idx = np.where(hgrd.mask_rho == 0)
#     h[idx] = hmin

#     # save raw bathymetry
#     hraw = h.copy()

#     # -- check bathymetry -- #
#     # check bathymetry roughness
#     RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
#     print('Max Roughness value is: ', RoughMat.max())

#     # # # smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
#     # rx0_max = 0.4
#     # h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

#     # check bathymetry roughness again
#     RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
#     print('Max Roughness value is: ', RoughMat.max())

#     return h,hraw


# def rotate_coords(xm, ym, ang_rot, degrees=False):
#     xm_mean = np.mean(xm)
#     ym_mean = np.mean(ym)

#     xm1 = xm.copy() - xm_mean
#     ym1 = ym.copy() - ym_mean

#     ang_rot_rad = ang_rot
#     if degrees:
#         ang_rot_rad = ang_rot/180 *np.pi

#     xrot = xm1*np.cos(ang_rot_rad) - ym1*np.sin(ang_rot_rad)
#     yrot = xm1*np.sin(ang_rot_rad) + ym1*np.cos(ang_rot_rad)

#     xrot += xm_mean
#     yrot += ym_mean
#     return xrot, yrot


# def interpolate_bathymetry(ddir, hgrd):
#     window = np.ones([5,5])
#     window = window/window.sum()
#     nc = xr.open_dataset(ddir)
#     topo = nc['elevation'].values #  np.loadtxt(os.path.join(ddir, 'etopo20data.gz'))
#     topo1 = signal.convolve2d(topo, window, 'same')
#     topo = topo1[::5,::5]
#     lons = nc['lon'][::5].values  # np.loadtxt(os.path.join(ddir, 'etopo20lons.gz'))
#     lats = nc['lat'][::5].values  # np.loadtxt(os.path.join(ddir, 'etopo20lats.gz'))

#     lonm, latm = np.meshgrid(lons,lats)
#     lonm = lonm.ravel()
#     latm = latm.ravel()

#     # depth positive
#     topo = -topo

#     # fix minimum depth
#     hmin = 5
#     topo = np.where(topo < hmin, hmin, topo)

#     # interpolate new bathymetry
#     lon, lat = np.meshgrid(lons, lats)
#     h = griddata((lonm.ravel(),latm.ravel()),topo.ravel(),(hgrd.lon,hgrd.lat), method='linear')
#     return h


# def quadrilateral_centroid(x,y):
#     # x = vertices[:, 0]
#     # y = vertices[:, 1]
#     x_shift = np.roll(x, 1)
#     y_shift = np.roll(y, 1)
#     area = 0.5*np.sum(x*y_shift - y*x_shift)
#     cx = (1/(6*area))*np.sum((x + x_shift)*(x*y_shift - y*x_shift))
#     cy = (1/(6*area))*np.sum((y + y_shift)*(x*y_shift - y*x_shift))
#     return np.array([cx, cy])


# def build_grid(x0,x1,y0,y1,lon0, lat0, rot, dxoffset,dyoffset):
#     # Define the projection
#     # projection = ccrs.AzimuthalEquidistant(central_latitude=lat0, central_longitude=lon0)

#     #  -- grid rotation -- #
#     # x = np.arange(-44,-40,1)+1
#     # original rectangle
#     x = np.array([x0, x1])# + xoffset
#     y = np.array([y0,y1])# + yoffset

#     # rectangle vertices
#     polygon = np.array([[x0,y0],[x0,y1],[x1,y1],[x1,y0]])

#     xm = polygon[:,0]
#     ym = polygon[:,1]

#     # rotation
#     xrot, yrot =rotate_coords(xm, ym, rot)

#     # convert to distance in meters
#     projection = pyproj.Proj(proj='aeqd', lon_0=lon0 + dxoffset,
#                                           lat_0=lat0 + dyoffset,
#                                           datum='WGS84')
#     x, y = projection(xrot, yrot)

#     # x, y = projection.transform_points(ccrs.PlateCarree(), xrot, yrot)[:, :2].T
#     centroid = quadrilateral_centroid(x,y)

#     # -- define diagonals of the quadrilater -- #
#     A = np.array([(x[1]+ x[2])/2, (y[1] + y[2])/2])
#     B = np.array([(x[3]+ x[0])/2, (y[3] + y[0])/2])
#     C = np.array([(x[0]+ x[1])/2, (y[0] + y[1])/2])
#     D = np.array([(x[2]+ x[3])/2, (y[2] + y[3])/2])

#     BA = B-A
#     DC = C-D

#     # get distance
#     distBA = (BA[0]**2 + BA[1]**2)**0.5
#     distDC = (DC[0]**2 + DC[1]**2)**0.5

#     # get size
#     nBA = int(distBA/3000)
#     nDC = int(distDC/3000)

#     # build a grid based on the diagonals
#     diag0 = np.linspace(-distBA/2,distBA/2, nBA) + centroid[0]
#     diag1 = np.linspace(-distDC/2,distDC/2, nDC) + centroid[1] 

#     # rotate the grid back and adjust using offsets
#     dxm, dym = rotate_coords(*np.meshgrid(diag0, diag1), -rot)
#     dxm = dxm #+ dxoffset
#     dym = dym #+ dyoffset

#     dxm = np.flipud(dxm).T
#     dym = np.flipud(dym).T
    

#     # Define the projection and rotate the coordinates
#     projection = pyproj.Proj(proj='aeqd', lon_0=lon0, lat_0=lat0, datum='WGS84')
#     # xrot, yrot = projection(x, y)

#     # Invert the transformation to get the original coordinates
#     xinvert, yinvert = pyproj.transform(projection, pyproj.Proj(proj='latlong', datum='WGS84'), dxm, dym)

#     # fig = plt.figure()
#     # ax = fig.add_subplot(projection=ccrs.PlateCarree())
#     # ax.coastlines()
#     # ax.scatter(xinvert,yinvert, marker='.', s=1)
#     # ax.plot(xinvert[0,0], yinvert[0,0],c='k', marker='*')
#     # ax.plot(xinvert[:,0], yinvert[:,0],c='r', label='xi_rho')
#     # ax.plot(xinvert[0,:], yinvert[0,:],c='b', label='eta_rho')
#     # plt.legend()
#     # Print the results
#     return xinvert, yinvert


# if __name__ == '__main__':
#     # plt.close('all')

#     # -- gets  the information from the config file -- #
#     reference = sys.argv[1] if len(sys.argv)==2  else 'swatl_2022_deep2'

#     dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
#     dicts = dicts[reference]
    
#     # -- paths -- #s
#     odir  = dicts['grid.output'] 
#     bfile = dicts['bathy_file']

#     # -- horizontal grid parameters -- #
#     dxdy        = dicts['grid.dxdy']
#     x0,x1,y0,y1 = dicts['grid.WESN']
#     xoffset     = dicts['grid.xoffset']
#     yoffset     = dicts['grid.yoffset']
#     rot         = dicts['grid.rot']

#     # -- vertical grid parameters -- #
#     N       = dicts['grid.N'] 
#     theta_s = dicts['grid.theta_s']
#     theta_b = dicts['grid.theta_b']
#     Tcline  = dicts['grid.Tcline']

#     # grid swatl_deep3
#     dxoffset=xoffset #* 1e5
#     dyoffset=yoffset #* 1e5
#     # rot = -50/180*np.pi
#     lon0 = -44.8
#     lat0 = -25.5


#     xrot, yrot = build_grid(x0,x1,y0,y1,lon0, lat0, rot, dxoffset,dyoffset)
#     # xrot = xrot.T
#     # yrot = yrot.T
#     hgrd, base = hgrid(xrot, yrot)

#     #pyroms.grid.edit_mask_mesh(hgrd, proj=base)
#     #edit_mask_mesh_ij is a faster version using imshow... but no base projection.
#     coast = pyroms.utility.get_coast_from_map(base)
#     # pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

#     topo = interpolate_bathymetry(bfile, hgrd)

#     # topo = -hraw1    # getting bathymetry (in this case the bathymetry was in the hromsg)
#     lon = xrot  # getting position
#     lat = yrot  # getting position
#     h,hraw = h_bathymetry(topo, lon, lat, hgrd)  # bathymetry [object (?)]

#     plt.close('all')
#     # vertical coordinate
#     vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
#     # print((vgrd.Cs_r*4000).astype('int'))
#     # plt.plot(vgrd.s_rho, vgrd.Cs_r, marker='x')


#     # ROMS grid
#     grd_name = 'BaciaSantos'
#     grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
#     fout = osp.join(odir)
#     # write grid to netcdf file
#     pyroms.grid.write_ROMS_grid(grd, filename=fout)
#     print(f'{fout} saved')


