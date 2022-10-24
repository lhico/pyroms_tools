import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.backend_bases import MouseButton
import xesmf as xe



def scaloa(xxr,yr,zmr, xg,yg, corrlenx, corrleny, err):
    # grid
    xi = xg.copy()
    yi = yg.copy()
    # phii = phig.copy()

    # observation
    x  = xxr.copy()
    y  = yr.copy()
    # phi  = phir.copy()
    t = zmr.copy()
    n=len(x)
    # squared distance matrix between observations
    scale =  (x[:,None] - x[None,:])**2/corrlenx**2
    scale +=  (y[:,None] - y[None,:])**2/corrleny**2
    A = (1 - err) * np.exp(-(scale))
    del(scale)
    # scale +=  abs(phi[:,None] - phi[None,:])/(phi[:,None]**2+phi[None,:]**2)**0.5 /corrlenphi


    # squared distance matrix between observation and grid
    if len(xi.shape) != 1:
        xc = xi.ravel()
        yc = yi.ravel()
        # phic = phii.ravel()
    else:
        xc = xi[:]
        yc = yi[:]
        # phic = phii[:]
    
    scalec =  (xc[:,None] - x[None,:])**2/corrlenx**2
    scalec +=  (yc[:,None] - y[None,:])**2/corrleny**2
    C = (1 - err) * np.exp(-(scalec))
    del scalec


    A = A + err * np.eye(len(A))

    tp = None
    ep = 1 - np.sum(C.T * np.linalg.solve(A, C.T), axis=0) / (1 - err)
    if any(t)==True: ##### was t!=None:
        t = t.ravel()
        tp = np.dot(C, np.linalg.solve(A, t))

    tp.reshape(xi.shape)
    return tp, ep



def adjust_bathymetry_between_grids(ncref, ncchange, corrlen, in_width, fr_width, err=0.01, varb='h', hmax=None):
    """
    Adjust bathymetry between grid by combining bathymetries and using objective analysis
    (optimum interpolation) to interpolate both grids. The idea is that the donor grid will
    be interpolated onto the external cells of the receiver grid and the bathymetry
    discrapencies will be smoothed and combined by the objective analysis.
    
    Args:
        ncref (xr.Dataset): Dataset with the donor bathymetry
        ncchange (xr.Dataset): Dataset with the receiver bathymetry
        corrlen (float): correlation length used to normalize the gaussian function (in degrees)
        in_width (integer): width of the external frame
        fr_width (integer): width of the external frame and a few cells of the receiver dataset
                            its value should be greater than `in_width`. This info will be used
                            in the bathymetry combination
        err (float, optional): Estimate of the normalized error. Defaults to 0.01.
        varb (str, optional): Name of the bathymetry dataarray. Defaults to 'h'.

    Returns:
        _type_: _description_
    """
    # we use xesmf as interpolation tool here 

    # 1) first rename dimension to xesmf standards
    # 2) interpolate from coarser to finer using xesmf
    # 3) replace 'interior' values of interpolated grid with ncchange values
    # 4) remove data in excess (objective analysis may need a LOT of memory)
    # 5) interpolate
    # 6) combine interpolated values and ncchange values

    # -- 1 -- #
    # renaming dimensions to xesmf standards
    ncchange = ncchange.rename(lon_rho='lon', lat_rho='lat')
    ncref = ncref.rename(lon_rho='lon', lat_rho='lat')
    ncchange = ncchange.set_coords(['lon', 'lat'])
    ncref = ncref.set_coords(['lon', 'lat'])

    # improving readability
    xmcha, ymcha = ncchange.lon.values, ncchange.lat.values

    # -- 2 -- #
    # interpolated coarse grid onto fine grid
    regridder = xe.Regridder(ncref, ncchange, 'nearest_s2d', extrap_method='nearest_s2d')
    interpolated = regridder(ncref)  # interpolating

    # -- 3 -- #
    # the zeroth and final 
    notborder = slice(in_width,-in_width)  # index slice where matrix will be nan

    # we are combining the interpolated and original bathymetries of the low resolution 
    # and high resolution grid into the highres grid
    # high resolution will be interior (notborder) points and low resolution will be
    # on the borders
    interpolated[varb].values[notborder,notborder] = ncchange.h[notborder,notborder] # set original values

    interpolated = updating_masks(interpolated, 0, uv=False)

    interpolated[varb].values[interpolated['mask_rho'].values==0] = np.nan
    # interpolated[varb].values[interpolated[varb].values==5] = np.nan   # erasing values considering depths

    # -- 4 -- #    
    # areas that are not the frame used in the bathymetry adjustment are removed
    interpolated[varb].values[fr_width:-fr_width,fr_width:-fr_width] = np.nan


    h = interpolated[varb].copy()
    # h[:,:3] = np.nan
    # h[:3,:] = np.nan
    # h[:,-3:] = np.nan
    # h[-3:,:] = np.nan
    isnotnan = ~np.isnan(h)
    

    # -- 5 -- #    
    # we adjust the bathymetry to the combined lowres and highres bathymetry present in the frame
    # using objective analysis 
    tp, err = scaloa(xmcha[isnotnan], ymcha[isnotnan],
                     h.values[isnotnan],
                     xmcha[isnotnan].ravel(), ymcha[isnotnan].ravel(),
                     corrlen, corrlen, err)

    interpolated = regridder(ncref)  # interpolating
    # -- 6 -- #    
    h1 = interpolated[varb].copy()
    # h1.values = tp.reshape(xmcha.shape)
    h1.values[isnotnan] = tp
    h1.values[fr_width:-fr_width, fr_width:-fr_width] = ncchange.h[fr_width:-fr_width, fr_width:-fr_width]
    
    if hmax is not None:
        h1.values[h1>hmax] = hmax

    return h1
    


def estimate_vertices(ncchange):
    """The netcdf created by coarse2fine.m does not include the vertices grid.
    We are adding it here. Note that this function is not general and depends on the
    grid rotation. 

    Args:
        ncchange (xr.Dataset): roms grid file (the function does not work as a subroutine)

    Returns:
        nc (xr.Dataset): roms grid file
    """
    nc = ncchange.copy()
    nc = nc.drop(['x_vert','y_vert','lon_vert','lat_vert'])
    print('WARNING, THIS FUNCITION IS NOT GENERAL. IT DEPENDS ON THE ROTATION OF THE GRID')
    shp = nc.x_psi.shape
    x_vert = np.zeros([shp[0]+2, shp[1]+2])
    x_vert[1:-1, 1:-1]  = nc.x_psi
    x_vert[0,:] = x_vert[1,:] - abs(x_vert[2,:]-x_vert[1,:])
    x_vert[-1,:] = x_vert[-2,:] + abs(x_vert[-3,:]-x_vert[-2,:])
    x_vert[:,0] = x_vert[:,1] - abs(x_vert[:,2]-x_vert[:,1])
    x_vert[:,-1] = x_vert[:,-2] + abs(x_vert[:,-2]-x_vert[:,-3])

    y_vert = np.zeros([shp[0]+2, shp[1]+2])
    y_vert[1:-1, 1:-1]  = nc.y_psi
    y_vert[0,:] = y_vert[1,:] - abs(y_vert[2,:]-y_vert[1,:])
    y_vert[-1,:] = y_vert[-2,:] + abs(y_vert[-3,:]-y_vert[-2,:])
    y_vert[:,0] = y_vert[:,1] + abs(y_vert[:,2]-y_vert[:,1])
    y_vert[:,-1] = y_vert[:,-2] - abs(y_vert[:,-2]-y_vert[:,-3])


    lon_vert = np.zeros([shp[0]+2, shp[1]+2])
    lon_vert[1:-1, 1:-1]  = nc.lon_psi
    lon_vert[0,:] = lon_vert[1,:] - abs(lon_vert[2,:]-lon_vert[1,:])
    lon_vert[-1,:] = lon_vert[-2,:] + abs(lon_vert[-3,:]-lon_vert[-2,:])
    lon_vert[:,0] = lon_vert[:,1] - abs(lon_vert[:,2]-lon_vert[:,1])
    lon_vert[:,-1] = lon_vert[:,-2] + abs(lon_vert[:,-2]-lon_vert[:,-3])

    lat_vert = np.zeros([shp[0]+2, shp[1]+2])
    lat_vert[1:-1, 1:-1]  = nc.lat_psi
    lat_vert[0,:] = lat_vert[1,:] - abs(lat_vert[2,:]-lat_vert[1,:])
    lat_vert[-1,:] = lat_vert[-2,:] + abs(lat_vert[-3,:]-lat_vert[-2,:])
    lat_vert[:,0] = lat_vert[:,1] + abs(lat_vert[:,2]-lat_vert[:,1])
    lat_vert[:,-1] = lat_vert[:,-2] - abs(lat_vert[:,-2]-lat_vert[:,-3])

    nc['lat_vert'] = (['eta_vert', 'xi_vert'], np.ma.MaskedArray(lat_vert,mask=False))
    nc['lon_vert'] = (['eta_vert', 'xi_vert'], np.ma.MaskedArray(lon_vert,mask=False))
    nc['x_vert'] = (['eta_vert', 'xi_vert'], np.ma.MaskedArray(x_vert,mask=False))
    nc['y_vert'] = (['eta_vert', 'xi_vert'], np.ma.MaskedArray(y_vert,mask=False))
    return nc


def updating_masks(ncchange, fr_width=0, hthreshold=5., uv=True):
    """
    Updates mask using the bathymetry as a reference.
    This scripts should be used after the update of bathymetry values

    Args:
        ncchange (xr.Dataset): ROMS grid
        fr_width (int): width of the frame where the hthreshold will be applied.
        hthreshold (float, optional): Depth threshold to set the landmask.
        uv (bool, optional): apply the update to u and v masks

    Returns:
        _type_: _description_
    """    
    nc = ncchange.copy()  # avoid changes of the original dataset
    nc.load()             # loading allows to update dataarray .values
    aux = (nc.h.values > hthreshold).astype(float)  # seamask 
    nc.mask_rho.values = aux 
    nc.mask_rho.values[fr_width:-fr_width] = nc.h.values[fr_width:-fr_width] > hthreshold
    # roms uses integer instead of booleans to deal with masks
    nc.mask_rho.values[fr_width:-fr_width] = nc.mask_rho.values[fr_width:-fr_width].astype(int)
    if uv:
        # updating u and v masks
        nc.mask_u.values = nc.mask_rho[:,1:]*nc.mask_rho[:,:-1]
        nc.mask_v.values = nc.mask_rho[1:,:]*nc.mask_rho[:-1,:]
    return nc

# use xesmf_env
if __name__ == '__main__':

    masklist = ['h', 'mask_u', 'mask_v']
    lonlist  = ['lon_rho', 'lon_u', 'lon_v']
    latlist  = ['lat_rho', 'lat_u', 'lat_v']

    ncchange = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools/data/treated/bacia_santos_nest_smooth.nc')
    ncref    = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools/data/treated/bacia_santos.nc')
    nc0      = xr.open_dataset('/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools/data/treated/bacia_santos.nc')


    # we need to combine info of the low res grid and the high res grid.
    # we use information from both grids

    h = adjust_bathymetry_between_grids(ncref, ncchange, 0.25, 8, 15, hmax=2500)


    ncchange = ncchange.rename(lon_rho='lon', lat_rho='lat')
    ncref = ncref.rename(lon_rho='lon', lat_rho='lat')
    ncchange = ncchange.set_coords(['lon', 'lat'])
    ncref = ncref.set_coords(['lon', 'lat'])

    h.values[np.isnan(h)] = 5
    ncchange['h'].values = h.values
    ncchange['hraw'].values =[h.values]

    # ncchange = updating_masks(ncchange,)
    # ncchange = ncchange.rename(lon='lon_rho', lat='lat_rho')


    # IMPORTANT
    # combine.m did not calculate _vert variables
    # this is a fix, that is NOT general and requires rectangular grid
    ncchange = estimate_vertices(ncchange)  # 

    # IMPORTANT,
    # combine.m script is producing nan at the following variables
    # here, we reinstate the variables back to the gridtop
    ncchange = ncchange.assign_coords(s_rho=nc0.s_rho.values)
    ncchange = ncchange.assign_coords(s_w=nc0.s_w.values)
    ncchange.Cs_r.values = nc0.Cs_r.values
    ncchange.hc.values   = nc0.hc.values
    ncchange.Cs_w.values = nc0.Cs_w.values

    ncchange = ncchange.rename(lon='lon_rho', lat='lat_rho')

    # unmatching masks at the borders may cause blow-up
    # we are updating masks only at points within a certain range
    # from the borders second argument in updating_masks
    ncchange = updating_masks(ncchange, fr_width=0)

    # plot the masks
    fig, ax =plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    ax = ax.ravel()
    ncref.h.plot.contour(x='lon', y='lat', levels=np.arange(0,2000,25), cmap=plt.cm.jet, extend='both', add_colorbar=False, ax=ax[0])
    ncchange.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(0,2000,25), cmap=plt.cm.jet, extend='both', add_colorbar=False, ax=ax[0])


    ncref.h.plot(x='lon', y='lat', levels=np.arange(0,2000,1), cmap=plt.cm.jet, extend='both', add_colorbar=False, ax=ax[1])
    ncchange.h.plot(x='lon_rho', y='lat_rho', levels=np.arange(0,2000,1), cmap=plt.cm.jet, extend='both', add_colorbar=False, ax=ax[1])
    # plt.scatter(xmcha, ymcha, marker='.', c='k', s=1)

    ncref.h.plot(x='lon', y='lat', levels=np.arange(-200,200,1), cmap=plt.cm.RdBu_r, extend='both', add_colorbar=False, ax=ax[2])
    ncchange.h.plot(x='lon_rho', y='lat_rho', levels=np.arange(-200,200,1), cmap=plt.cm.RdBu_r, extend='both', add_colorbar=False, ax=ax[2])
    ncchange.h.plot.contour(x='lon_rho', y='lat_rho', levels=np.arange(0,2000,25), cmap=plt.cm.RdBu_r, extend='both', add_colorbar=False, ax=ax[3])


    ncchange.to_netcdf('/home/otel/Dropbox/trabalho_irado/2022/postdoc/202203_ceresIV/pyroms_tools/data/treated/bacia_santos_nest_adjusted.nc')


