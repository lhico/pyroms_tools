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
    parser.add_argument('--bathy_file', type=str, help='Bathymetry file path.', required=False)
    parser.add_argument('--dxdy', type=float, help='Grid spacing.', required=False)
    parser.add_argument('--wesn', type=float, nargs=4, metavar=('x0', 'x1', 'y0', 'y1'), help='Grid boundaries (WESN).', required=False)
    parser.add_argument('--xoffset', type=float, help='X offset.', required=False)
    parser.add_argument('--yoffset', type=float, help='Y offset.', required=False)
    parser.add_argument('--rot', type=float, help='Grid rotation angle.', required=False)
    parser.add_argument('--N', type=int, help='Number of vertical levels.', required=False)
    parser.add_argument('--theta_s', type=float, help='Surface control parameter.', required=False)
    parser.add_argument('--theta_b', type=float, help='Bottom control parameter.', required=False)
    parser.add_argument('--Tcline', type=float, help='Thermocline depth.', required=False)
    parser.add_argument('--gridout', type=str, help='Output grid file path.', required=False)
    parser.add_argument('--reference', type=str, default='default', help='Reference key for the configuration.', required=False)
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
            # Process xa and ya as needed
        elif isinstance(geom, MultiPolygon):
            for poly in geom.geoms:  # Use geom.geoms to iterate over polygons in MultiPolygon
                coords = poly.exterior.coords.xy
                xa = np.array(coords[0], np.float32)
                ya = np.array(coords[1], np.float32)
                # Process xa and ya as needed

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
    
    if args.config:
        dicts = load_config_from_yaml(args.config, args.reference)
    else:
        dicts = {
            'bathy_file': args.bathy_file,
            'grid.dxdy': args.dxdy,
            'grid.WESN': args.wesn,
            'grid.xoffset': args.xoffset,
            'grid.yoffset': args.yoffset,
            'grid.rot': args.rot,
            'grid.N': args.N,
            'grid.theta_s': args.theta_s,
            'grid.theta_b': args.theta_b,
            'grid.Tcline': args.Tcline,
            'grid.grid': args.gridout
        }
    
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