import ast
import xarray as xr
import numpy as np
from scipy import interpolate as itp
import netCDF4 as nc
try:
    import pyroms
except:
    pass
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


class nudgcoef():
    ''' A class to write the Nudging coeficient file for ROMS '''

    def __init__(self,roms_grid):
        ''' init an object of the class with the pyroms grid ID ''' 
        self.grd = pyroms.grid.get_ROMS_grid(roms_grid)
        return None

    def __call__(self,east_dict,west_dict,north_dict,south_dict,tracer_timescales,foutname='./nudging_coef.nc'):
        ''' call with following dictionaries :
        4 boundaries dict + tracer timescales
        for example :
        east_dict  = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        west_dict  = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        north_dict = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        south_dict = {'nudge':True,'factor': 1,'width':50,'transition':'linear'}
        tracer_timescales = {'M2':30,'M3':30,'temp':30,'salt':30,'tracer':30}
        tips:
        * nudge = True if open boundary, False otherwise
        * factor allows to have different timescales at each boundary
        * width is in grid points
        * transition shapes how timescale varies spatially
        * tracer timescales are in days
        '''
        self.east_dict = east_dict
        self.west_dict = west_dict
        self.north_dict = north_dict
        self.south_dict = south_dict
        self.tra_ts = tracer_timescales
        self.foutname = foutname
        # create 2d coef
        self.nud2 = self._create_nudgcoef_2d()
        # create 3d coef
        self.nud3 = self._create_nudgcoef_3d()
        # write to netcdf
        self._write_nc_file()
        return None

    def _create_nudgcoef_3d(self):
        ''' expand 2d coef along the vertical '''
        # RD: later we could imagine multiplying by
        # a vertical profile if needed
        ny, nx = self.grd.hgrid.mask_rho.shape
        nz = self.grd.vgrid.N
        nudgcoef = np.zeros((nz,ny,nx))
        for kz in np.arange(nz):
            nudgcoef[kz,:,:] = self.nud2[:,:]
        return nudgcoef

    def _create_nudgcoef_2d(self):
        ''' create the 2d nudging coef from dictionaries '''
        ny, nx = self.grd.hgrid.mask_rho.shape
        nudgcoef_west  = np.zeros((ny,nx))
        nudgcoef_east  = np.zeros((ny,nx))
        nudgcoef_north = np.zeros((ny,nx))
        nudgcoef_south = np.zeros((ny,nx))
        nudgcoef       = np.zeros((ny,nx))
        mask           = self.grd.hgrid.mask_rho
        # west boundary
        if self.west_dict['nudge'] is True:
            fc = self.west_dict['factor']
            wd = self.west_dict['width']
            tr = self.west_dict['transition']
            if tr == 'linear':
                for ji in np.arange(0,wd):
                    nudgcoef_west[:,ji] = fc * (wd-ji) / float(wd)
            elif tr == 'linear_nocoast':
                for ji in np.arange(0,wd):
                    nudgcoef_west[:,ji] = mask[:,0] * fc * (wd-ji) / float(wd)
            else:
                print('transition not coded') ; pass
        # east boundary
        if self.east_dict['nudge'] is True:
            fc = self.east_dict['factor']
            wd = self.east_dict['width']
            tr = self.east_dict['transition']
            if tr == 'linear':
                for ji in np.arange(nx-wd,nx):
                    nudgcoef_east[:,ji] = fc * (wd-nx+ji) / float(wd)
            elif tr == 'linear_nocoast':
                for ji in np.arange(nx-wd,nx):
                    nudgcoef_east[:,ji] = mask[:,-1] * fc * (wd-nx+ji) / float(wd)
            else:
                print('transition not coded') ; pass
        # south boundary
        if self.south_dict['nudge'] is True:
            fc = self.south_dict['factor']
            wd = self.south_dict['width']
            tr = self.south_dict['transition']
            if tr == 'linear':
                for jj in np.arange(0,wd):
                    nudgcoef_south[jj,:] = fc * (wd-jj) / float(wd)
            if tr == 'linear_nocoast':
                for jj in np.arange(0,wd):
                    nudgcoef_south[jj,:] = mask[0,:] * fc * (wd-jj) / float(wd)
            else:
                print('transition not coded') ; pass
        # north boundary
        if self.north_dict['nudge'] is True:
            fc = self.north_dict['factor']
            wd = self.north_dict['width']
            tr = self.north_dict['transition']
            if tr == 'linear':
                for jj in np.arange(ny-wd,ny):
                    nudgcoef_south[jj,:] = fc * (wd-ny+jj) / float(wd)
            if tr == 'linear_nocoast':
                for jj in np.arange(ny-wd,ny):
                    nudgcoef_south[jj,:] = mask[-1,:] * fc * (wd-ny+jj) / float(wd)
            else:
                print('transition not coded') ; pass


        # create the total coefficient by combining all 4 fields
        # the max functions is useful to make nice corners when
        # individual field overlap
        # maybe not the most efficient but short and readable
        for jj in np.arange(ny):
            for ji in np.arange(nx):
                nudgcoef[jj,ji] = max(nudgcoef_west[jj,ji], \
                nudgcoef_east[jj,ji],nudgcoef_north[jj,ji],nudgcoef_south[jj,ji])
        return nudgcoef

    def _write_nc_file(self):
        ''' writing to netcdf and multiplying by inverse timescales '''
        ncfile = self.foutname
        fid = nc.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')

            # dimensions
        fid.createDimension('xi_rho', np.size(self.grd.hgrid.mask_rho,1))
        fid.createDimension('eta_rho', np.size(self.grd.hgrid.mask_rho,0))
        fid.createDimension('s_rho', self.grd.vgrid.N)
        fid.createDimension('s_w', self.grd.vgrid.Np)
        fid.description = 'Nudging coefficients for grid' + self.grd.name
        # vertical coordinate
        fid.createVariable('s_rho', 'f8', ('s_rho'))
        fid.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
        fid.variables['s_rho'].valid_min = '-1'
        fid.variables['s_rho'].valid_max = '0'
        fid.variables['s_rho'].field = 's_rho,scalar'
        fid.variables['s_rho'][:] = self.grd.vgrid.s_rho

            # variables
        O_M2_NudgeCoef     = fid.createVariable('M2_NudgeCoef',     'f8', ('eta_rho','xi_rho',))
        O_M3_NudgeCoef     = fid.createVariable('M3_NudgeCoef',     'f8', ('s_rho','eta_rho','xi_rho',))
        O_temp_NudgeCoef   = fid.createVariable('temp_NudgeCoef',   'f8', ('s_rho','eta_rho','xi_rho',))
        O_salt_NudgeCoef   = fid.createVariable('salt_NudgeCoef',   'f8', ('s_rho','eta_rho','xi_rho',))
        O_tracer_NudgeCoef = fid.createVariable('tracer_NudgeCoef', 'f8', ('s_rho','eta_rho','xi_rho',))
        # data
        O_M2_NudgeCoef[:,:]       = (1./self.tra_ts['M2'])     * self.nud2
        O_M3_NudgeCoef[:,:,:]     = (1./self.tra_ts['M3'])     * self.nud3
        O_temp_NudgeCoef[:,:,:]   = (1./self.tra_ts['temp'])   * self.nud3
        O_salt_NudgeCoef[:,:,:]   = (1./self.tra_ts['salt'])   * self.nud3
        O_tracer_NudgeCoef[:,:,:] = (1./self.tra_ts['tracer']) * self.nud3
        # attributes
        O_M2_NudgeCoef.long_name = '2D momentum inverse nudging coefficients'
        O_M2_NudgeCoef.units = 'days-1'
        O_M2_NudgeCoef.coordinates = 'xi_rho eta_rho'

        O_M3_NudgeCoef.long_name = '3D momentum inverse nudging coefficients'
        O_M3_NudgeCoef.units = 'days-1'
        O_M3_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_temp_NudgeCoef.long_name = 'temp inverse nudging coefficients'
        O_temp_NudgeCoef.units = 'days-1'
        O_temp_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_salt_NudgeCoef.long_name = 'salt inverse nudging coefficients'
        O_salt_NudgeCoef.units = 'days-1'
        O_salt_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'

        O_tracer_NudgeCoef.long_name = 'generic tracer inverse nudging coefficients'
        O_tracer_NudgeCoef.units = 'days-1'
        O_tracer_NudgeCoef.coordinates = 'xi_rho eta_rho s_rho'
        # close
        fid.close()
        return None