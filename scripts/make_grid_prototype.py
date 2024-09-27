import matplotlib.pyplot as plt
import os.path as osp
from scipy.interpolate import griddata
import xarray as xr
# from utils import configs
from scipy import ndimage
from scipy import signal
from utils import utils as ut
from dry_package.plot_schemes import maps
import numpy as np


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

if __name__ == '__main__':

    # grid_config_pyroms.txt is a python dictionary with the information
    # used in this script to prepare the grid

    # -- gets  the information from the config file -- #
    reference = 'pbs_202109_glorys'
    dicts = ut._get_dict_paths('../configs/grid_config_pyroms.txt')
    dicts = dicts[reference]

    # -- paths -- #s
    odir  = dicts['output_dir'] # directoy where the file will be saved
    
    # -- horizontal grid parameters -- #
    dxdy    = dicts['dxdy']
    x0,x1,y0,y1 = dicts['WESN']
    xoffset = dicts['xoffset']
    yoffset = dicts['yoffset']
    rot     = dicts['rot']

    # -- vertical grid parameters -- #
    N       = dicts['N'] 
    theta_s = dicts['theta_s']
    theta_b = dicts['theta_b']
    Tcline  = dicts['Tcline']


    #  -- grid rotation -- #
    # x = np.arange(-44,-40,1)+1
    x = np.arange(x0,x1,dxdy) + xoffset
    y = np.arange(y0,y1,dxdy) + yoffset


    xm, ym = np.meshgrid(x,y)
    xrot, yrot =rotate_coords(xm, ym, rot)
    # ax.scatter(xrot,yrot,s=1,marker='.')
    # -- define horizontal roms grid object -- #

    # plt.close('all')
    ax = maps.make_overall(extent=[-55,-30,-35,-15])
    ax.scatter(xrot,yrot, marker='.', s=1)
    ax.coastlines(resolution='10m')