import xarray as xr
import numpy as np
from tqdm import tqdm

###############################################################################


def conversion_helper():

    d = {
        'default_conversion': default_conversion,
        'Pa_to_mb': Pa_to_mb,
        'rate_tp_to_inst': rate_tp_to_inst,
        'rate_e_to_inst': rate_e_to_inst
    }

    return d

###############################################################################]


def Pa_to_mb(ds, varname='msl'):

    return ds[varname]/100

###############################################################################


def default_conversion(ds, varname, dt=1):
    """
    dt: hours
    varname: variable to convert
    ds: xarray.Dataset
    """
    return ds[varname]/(3600*dt)

###############################################################################


def rate_tp_to_inst(ds, varname, rho_w=1000, dt=1):

    return (ds[varname] * rho_w)/(dt*3600)

###############################################################################


def rate_e_to_inst(ds, varname, rho_w=1000, dt=1):
    """
    
    """
    return (ds[varname] * rho_w)/(dt*3600)

###############################################################################


def vapor_pressure(T):
    """
    e = vapor_pressure(T)

    actual vapor pressure in hPa (millibars) from dewpoint temperature `T` in degree Celsius
    using using [1]. If `T` is the air temperature, then  `e` is the saturated vapor
    pressure over liquid water is given by:

    ``
    e(T) = 6.11 \\cdot 10 ^ {\\left(  \\frac{7.5 T}{237.7 + T} \\right)}
    ``

    [1] https://web.archive.org/web/20200926200733/https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
    """
    return 6.11 * 10.0 ** (7.5 * T / (237.7 + T))

###############################################################################


def relative_humidity(t2m_C, d2m_C):
    """
    rh = relative_humidity(t2m_C,d2m_C)

    Compute the relative humidity (between 0 and 100) from temperature at 2 m, and dew_temperature at
    2 m) both in degree Celsius)

    [1] https://web.archive.org/web/20200926200733/https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
    """

    return 100 * vapor_pressure(d2m_C) / vapor_pressure(t2m_C)

###############################################################################


def add_global_attributes(ds):

    ds.attrs['history'] = 'Forcing file created with [script_name]'
    ds.attrs['type'] = 'Forcing file'
    ds.attrs['title'] = 'ECMWF ERA5 Dataset for South Brazil Bight (SBB)'

    return ds

###############################################################################

class dQdT(object):
    """ 
    Based on Barnier, Siefridt and Marchesiello (1995) available at:

    >> original citation: https://doi.org/10.1016/0924-7963(94)00034-9
    >> open access at >> https://sci-hub.se/https://doi.org/10.1016/0924-7963(94)00034-9
    """
    def __init__(self, scalarwind, pressure, T, rhoair=1):
        self.T = T              # sea surface temperature [K]
        self.Ua = scalarwind    # scalar wind speed at anemometer level [m s-1]
        self.rhoair = rhoair    # air density [kg m-3]
        self.Pa = pressure      # mean sea level pressure [mb]


    def dQdT_infrared(self):
        """Equation 5a"""
        
        sigma = 5.67 * 10**-8  # Steffan-Boltzmann constant  [W m**-2*K**-4]
        self.dQdT_ir = -4 * sigma * self.T**3


    def dQdT_sensibleheat(self):
        """Equation 5b"""
        Cp = 1.0048 * 10**3  # air specific heat [J kg**-1 **K-1 ]
        Ch = 10**-3          # Bulk transfer coeficient for sensible heat
        self.dQdT_sh = - self.rhoair * Cp * Ch * self.Ua


    def dQdT_latentheat(self):
        """Equation 5c"""

        Ce = 1.15e-3
        L = 2.508 * 1e6                         # latent heat of vaporization [J kg**-1]
        es = lambda T: 10 ** (9.4051 - 2353/T)  # saturated water pressure vapo
        qs = lambda p,T: 0.622/p * es(T)        # humidity

        self.dQdT_lh = -self.rhoair * Ce * L * self.Ua * \
                       2353*np.log(10)*qs(self.Pa, self.T)/self.T**2

    def dQdT(self):
        """See equation 6b"""
        self.dQdT_infrared()
        self.dQdT_latentheat()
        self.dQdT_sensibleheat()
        return self.dQdT_lh + self.dQdT_ir + self.dQdT_sh

###############################################################################

def scaling(x):
    sf = x.attrs['scale_factor']
    ao = x.attrs['add_offset']
    return sf*x + ao

###############################################################################

# internal function to extrapolate using xESMF function
def _extrapolate_xesmf(dsout, dst, mask, varname, method_interp='bilinear', method_extrap='inverse_dist'):
    """
    """
    import xesmf as xe

    dsout['mask'] = mask
    regridder = xe.Regridder(dsout, dst, method=method_interp, extrap_method=method_extrap)
    dsout = regridder(dsout[varname]).to_dataset(name=varname)

    return dsout

# internal function to extraploate using the Laplace operator
def _extrapolate_laplace(src, mask, varb):
    """
    """
    import utils.extrapolate as ex

    # variable to extrapolate
    varz = src[varb].values

    # applying Laplace's method to fill nearshore gaps
    undef = 2.0e+35
    tx = 0.9 * undef
    critx = 0.01
    cor = 1.6
    mxs = 10

    nx = src.dims['longitude']
    ny = src.dims['latitude']
    toxi = nx
    toeta = ny

    # looping through time to apply the extrapolation
    field = np.copy(varz)
    for t in tqdm(range(field.shape[0])):
        # extracting 2D field to extrapolate
        tmp_field = field[t,:,:]

        # masking based on SST landmask
        tmp_field = np.ma.masked_where(mask, tmp_field)

        # replacing mask by undef value so the extrapolator can work over it
        tmp_field = np.where(np.ma.getmask(tmp_field), undef, tmp_field)


        field[t,:,:] = ex.extrapolate.fill(int(1),int(toxi),
                                           int(1),int(toeta),
                                           float(tx), float(critx), float(cor), float(mxs),
                                           np.asarray(tmp_field, order='F'),
                                           int(nx),
                                           int(ny))

    dsout = src[[varb]].copy(deep=True)
    dsout[varb] = (("time", "latitude", "longitude"), field)

    return dsout

# main function to extrapolate the data
def extrapolating_era5(src, varb, sst, extrapolate_method='laplace', dst=None, mask=None):
    """
    This function only extrapolate ocean data onto the land, using the SST landmask from the ERA5 dataset. It is important
    to mention that this is not an interpolation to the ROMS grid, therefore, the ERA5 data covers, and must do, a much broader 
    area than the domain modeled, so ROMS can internally interpolate.
    """
    # creating landmask based on SST
    #mask = xr.where(src['sst'].isel(time=0) == src['sst'].attrs['missing_value'], 1, 0).values
    if mask is None:
        # if mask is not sent, then create based on ERA5 SST. Note that to use the extrapolation
        # method from xesmf, you must send a mask compatible with the destiny grid
        mask = xr.where(sst.isel(time=0) == sst.attrs['missing_value'], 1, 0).values

    if extrapolate_method == 'laplace':
        dsout = _extrapolate_laplace(src, mask, varb)
    elif extrapolate_method == 'xesmf' and dst:
        # using the same domain as sourc and destiny, because we want to keep the ERA5 grid so ROMS will interpolate internally to its domain
        dsout = _extrapolate_xesmf(src, src, mask, varb)
    elif extrapolate_method == 'xesmf' and not dst:
        raise ValueError('You must provide a destiny grid (dst) in order to use the xESMF extrapolation method')
    else:
        raise ValueError(f"Method unavailable ({extrapolate_method}. Please select a valid option: laplaca or xesmf")

    return dsout
