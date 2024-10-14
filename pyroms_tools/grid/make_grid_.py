import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import xarray as xr
import pyroms
from scipy import signal
from bathy_smoother import bathy_tools, bathy_smoothing
import yaml
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import logging
from pyproj import Proj
from shapely.geometry import Polygon, MultiPolygon


# Constants
HMIN = 5

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate grid configuration for pyroms.')
    parser.add_argument('--config', type=str, help='Path to the YAML configuration file.')

    return parser.parse_args()

def load_config_from_yaml(yaml_file, reference):
    """Load configuration from a YAML file."""
    try:
        with open(yaml_file, 'r') as file:
            config = yaml.safe_load(file)
        return config[reference]
    except Exception as e:
        logging.error(f"Error loading YAML file: {e}")
        sys.exit(1)

def hgrid(lon_rho, lat_rho):
    """Create horizontal grid and mask land areas."""
    lonv = lon_rho.copy()
    latv = lat_rho.copy()
    # projection = ccrs.Mercator()
    projection = Proj(proj='merc', ellps='WGS84')
    hgrd = pyroms.grid.CGrid_geo(lonv, latv, projection)
    
    # Create a mask using Cartopy's Natural Earth features
    land = cfeature.NaturalEarthFeature('physical', 'land', '10m')
    for geom in land.geometries():
        if isinstance(geom, Polygon):
            coords = geom.exterior.coords.xy
            xa = np.array(coords[0], np.float32)
            ya = np.array(coords[1], np.float32)
            xa, ya = hgrd.proj(xa, ya)
            vv = np.zeros((xa.shape[0],2))
            vv[:, 0] = xa
            vv[:, 1] = ya
            hgrd.mask_polygon(vv,mask_value=0)
            # Process xa and ya as needed
        elif isinstance(geom, MultiPolygon):
            for poly in geom.geoms:  # Use geom.geoms to iterate over polygons in MultiPolygon
                coords = poly.exterior.coords.xy
                xa = np.array(coords[0], np.float32)
                ya = np.array(coords[1], np.float32)
                xa, ya = hgrd.proj(xa, ya)

                vv = np.zeros((xa.shape[0],2))
                vv[:, 0] = xa
                vv[:, 1] = ya
                hgrd.mask_polygon(vv,mask_value=0)

    return hgrd, projection

def h_bathymetry(topo, lon, lat, hgrd):
    """Interpolate bathymetry data onto the grid."""
    topo = np.where(topo < HMIN, HMIN, topo)
    h = griddata((lon.ravel(), lat.ravel()), topo.ravel(), (hgrd.lon_rho, hgrd.lat_rho), method='linear')
    h = np.where(h < HMIN, HMIN, h)
    idx = np.where(hgrd.mask_rho == 0)
    h[idx] = HMIN
    hraw = h.copy()
    RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
    logging.info(f'Max Roughness value is: {RoughMat.max()}')
    return h, hraw

def rotate_coords(xm, ym, ang_rot, degrees=False):
    """Rotate coordinates by a given angle."""
    xm_mean = np.mean(xm)
    ym_mean = np.mean(ym)
    xm1 = xm.copy() - xm_mean
    ym1 = ym.copy() - ym_mean
    ang_rot_rad = ang_rot if not degrees else ang_rot / 180 * np.pi
    xrot = xm1 * np.cos(ang_rot_rad) - ym1 * np.sin(ang_rot_rad)
    yrot = xm1 * np.sin(ang_rot_rad) + ym1 * np.cos(ang_rot_rad)
    xrot += xm_mean
    yrot += ym_mean
    return xrot, yrot

def interpolate_bathymetry(ddir, hgrd):
    """Interpolate bathymetry data from a file."""
    window = np.ones([5, 5]) / 25
    nc = xr.open_dataset(ddir)
    topo = nc['elevation'].values
    topo = signal.convolve2d(topo, window, 'same')[::5, ::5]
    lons = nc['lon'][::5].values
    lats = nc['lat'][::5].values
    lonm, latm = np.meshgrid(lons, lats)
    lonm = lonm.ravel()
    latm = latm.ravel()
    topo = -topo
    topo = np.where(topo < HMIN, HMIN, topo)
    lon, lat = np.meshgrid(lons, lats)
    h = griddata((lonm.ravel(), latm.ravel()), topo.ravel(), (hgrd.lon, hgrd.lat), method='linear')
    return h

def main():
    """Main function to generate the grid configuration."""
    args = parse_args()
    
    dicts = load_config_from_yaml(args.config, 'default')

    
    # Extract parameters from the configuration
    bfile = dicts['bathy_file']
    dxdy = dicts['grid']['dxdy']
    x0, x1, y0, y1 = dicts['grid']['WESN']
    xoffset = dicts['grid']['xoffset']
    yoffset = dicts['grid']['yoffset']
    rot = dicts['grid']['rot']
    N = dicts['grid']['N']
    theta_s = dicts['grid']['theta_s']
    theta_b = dicts['grid']['theta_b']
    Tcline = dicts['grid']['Tcline']
    gridout = dicts['grid']['grid']

    # Generate grid coordinates
    x = np.arange(x0, x1, dxdy) + xoffset
    y = np.arange(y0, y1, dxdy) + yoffset
    xm, ym = np.meshgrid(x, y)
    xrot, yrot = rotate_coords(xm, ym, rot)
    
    # Create horizontal grid and interpolate bathymetry
    hgrd, projection = hgrid(xrot, yrot)
    topo = interpolate_bathymetry(bfile, hgrd)
    lon = xrot
    lat = yrot
    h, hraw = h_bathymetry(topo, lon, lat, hgrd)
    
    # Close all plots
    plt.close('all')
    
    # Create vertical grid and ROMS grid
    vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
    grd_name = 'BaciaSantos'
    grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)
    
    # Write the grid to a file
    pyroms.grid.write_ROMS_grid(grd, filename=gridout)

if __name__ == '__main__':
    main()