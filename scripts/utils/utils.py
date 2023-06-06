import ast
import xarray as xr
import numpy as np
from scipy import interpolate as itp
# from metpy import interpolate as metitp


def _get_dict_paths(fpath):
    """Read a text file with a dictionary and returns
    the dictionary"""
    file = open(fpath, 'r')
    contents = file.read()
    dictionary = ast.literal_eval(contents)
    return dictionary

def get_ref_dicts(fpath, ref, src_path, out_path):
    """Read and select information from the config file (fpath)
    and set the appropriate paths into a python dictionary.

    The config file contains dictionaries with useful paths, boundary
    limits, grid resolution, among other info.

    Args:
        fpath (string): path to the config file
        ref (string)  : this key selects a dataset from the config file
        src_path (string): this string will replace the <srcpath> in
                           appropriate items from the config file.
        out_path (string): this string will replace the <outpath> in
                           appropriate items from the config file.
             
    Returns:
        [dict]: dictionary with the selected item
    """    
    dicts = _get_dict_paths(fpath)
    dicts = dicts[ref]
    temp = {}
    for bla in dicts.keys():
        if type(dicts[bla]) is str:
            temp[bla] = dicts[bla].replace('<srcpath>', src_path)
            temp[bla] = temp[bla].replace('<outpath>', out_path)
        else:
            temp[bla] = dicts[bla]
    print(temp)
    return temp

# def _nn_interpolation(xc,yc,zc,xo,yo):
#     """Wrapper for metpy natural neighbor interpolation method.

#     Args:
#         xc (np.ndarray): x input grid component(2D array)
#         yc (np.ndarray): y input grid component(2D array)
#         zc (np.ndarray): grodded datat(2D array)
#         xo (np.ndarray): x output grid component(2D array)
#         yo (np.ndarray): y output grid component(2D array)
#         offset (float, optional): metpy has an issue with natural neighbor
#                 gridding, when the interpolation points coincide with the source
#                 grid. We add this offset as a quick-fix. Defaults to 1e-6. See:
#                 https://github.com/Unidata/MetPy/issues/346

#     Returns:
#         [type]: [description]
#     """    

#     xyp = np.array([xc.ravel(), yc.ravel()]).T
#     xyo = np.array([xo.ravel(), yo.ravel()]).T
    
#     zo = metitp.natural_neighbor_to_points(xyp, zc.ravel(), xyo)

#     return zo.reshape(xo.shape)


def _nn_interpolation(xc,yc,zc,xo,yo):
    """Wrapper for scipy.interpolate.griddata

    Args:
        xc (np.ndarray): x input grid component(2D array)
        yc (np.ndarray): y input grid component(2D array)
        zc (np.ndarray): gridded data(2D array)
        xo (np.ndarray): x output grid component(2D array)
        yo (np.ndarray): y output grid component(2D array)
    Returns:
        [type]: [description]
    """    
    # we 
    xyp = np.array([xc.ravel(), yc.ravel()]).T
    xyo = np.array([xo.ravel(), yo.ravel()]).T
    
    zo = itp.griddata(xyp, zc.ravel(), xyo, method='linear')

    return zo.reshape(xo.shape)


def refine_grid(nc, Gfactor=3):
    lon_rho = nc.lon_rho.values
    lat_rho = nc.lat_rho.values
    # lon_psi = nc.lon_psi.values
    # lat_psi = nc.lat_psi.values
    # lon_u = nc.lon_u.values
    # lat_u = nc.lat_u.values
    # lon_v = nc.lon_v.values
    # lat_v = nc.lat_v.values

    # x_rho = nc.x_rho.values
    # y_rho = nc.y_rho.values
    # x_psi = nc.x_psi.values
    # y_psi = nc.y_psi.values
    # x_u = nc.x_u.values
    # y_u = nc.y_u.values
    # x_v = nc.x_v.values
    # y_v = nc.y_v.values

    # Gfactor = 3   # grid refinement factor
    
    [Lp, Mp] = np.shape(lon_rho)
    
    L = Lp-1
    M = Mp-1

    # set coarser grid fractional coordinates
    half = np.floor(Gfactor).astype(int)

    Imin = half     # Coarse grid lower-left  I-coordinate (PSI-point)
    Imax = Lp-half  # Coarse grid upper-right I-coordinate (PSI-point) 
    Jmin = half     # Coarse grid lower-left  J-coordinate (PSI-point) 
    Jmax = Mp-half  # Coarse grid upper-right J-coordinate (PSI-point) 

    Imin = gridfocus1[0][0]
    Imax = gridfocus1[0][1]
    Jmin = gridfocus1[1][0]
    Jmax = gridfocus1[1][1]


    YrC, XrC = np.meshgrid(np.arange(0.5,Mp+0.5)  , np.arange(0.5,Lp+0.5))
    # YpC, XpC = np.meshgrid(np.arange(1.0,M+1)     , np.arange(1,L+1))
    # YuC, XuC = np.meshgrid(np.arange(0.5,Mp-0.5)  , np.arange(1,L+2))
    # YvC, XvC = np.meshgrid(np.arange(1.0,M+2)     , np.arange(0.5,Lp-0.5))

    # set finer grid fractional coordinates
    delta     = 1.0 / Gfactor
    halfdelta = 0.5 * delta

    IpF = np.arange(Imin,Imax + delta, delta)
    JpF = np.arange(Jmin,Jmax + delta, delta)
    IrF = np.insert(IpF+halfdelta, 0, IpF[0]-halfdelta)
    JrF = np.insert(JpF+halfdelta, 0, JpF[0]-halfdelta)

    YrF, XrF = np.meshgrid(JrF, IrF)
    # YpF, XpF = np.meshgrid(JpF, IpF)
    # YuF, XuF = np.meshgrid(JrF, IpF)
    # YvF, XvF = np.meshgrid(JpF, IrF)

    xyp = np.array([XrC.ravel(), YrC.ravel()]).T
    # xyo = np.array([XrF.ravel(), YrF.ravel()]).T
    
    # interpolation
    lonrF = _nn_interpolation(XrC, YrC, lon_rho, XrF, YrF)
    latrF = _nn_interpolation(XrC, YrC, lat_rho, XrF, YrF)

    # lonpF = _nn_interpolation(XpC, YpC, lon_psi, XpF, YpF)
    # latpF = _nn_interpolation(XpC, YpC, lat_psi, XpF, YpF)

    # lonuF = _nn_interpolation(XuC, YuC, lon_u, XuF, YuF)
    # latuF = _nn_interpolation(XuC, YuC, lat_u, XuF, YuF)

    # lonvF = _nn_interpolation(XvC, YvC, lon_v, XvF, YvF)
    # latvF = _nn_interpolation(XvC, YvC, lat_v, XvF, YvF)

    return lonrF, latrF

def create_grid_file(id_val, name_val, grdfile_val, N_val, grdtype_val, Vtrans_val, theta_s_val, theta_b_val, Tcline_val):
    content = f'''#
id      = {id_val}
name    = {name_val}
grdfile = {grdfile_val}
N       = {N_val}
grdtype = {grdtype_val}
Vtrans  = {Vtrans_val}
theta_s = {theta_s_val}
theta_b = {theta_b_val}
Tcline  = {Tcline_val}
#'''

    with open('gridid.txt', 'w') as file:
        file.write(content)


def update_mask(ds):
    ds = ds.copy()
    ds['mask_u'].values = ds['mask_rho'][:,:-1].values * ds['mask_rho'][:,1:].values
    ds['mask_v'].values = ds['mask_rho'][:-1,:].values * ds['mask_rho'][1:,:].values
    ds['mask_psi'].values = ds.mask_u[:-1, :].values * \
                        ds.mask_v[:, :-1].values
    return ds
