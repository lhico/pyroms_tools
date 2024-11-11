#!/usr/bin/python2.7
# -*-coding:Utf-8-*
import numpy as np
from matplotlib.path import Path
import _pickle as pickle
import datetime
import scipy.interpolate as interpolate
import os.path
import xesmf
from netCDF4 import Dataset

#__author__   = 'Thomas Roc'
#__email__    = 'thomas.roc@gmail.co.uk'
#__contributors__ = 'Andressa Palma, Danilo Silva, Dalton Sasaki'
#__email__    = 'thomas.roc@gmail.co.uk'
#__created__  = datetime.datetime(2014, 12, 5)
#__modified__ = datetime.datetime(2023, 09, 17)
#__version__  = "0.0"
#__status__   = "Under Development"

def read_tpxo_data(lon,lat, harmonic_names,tpxofile_data,
                   pickle_name='tpxo8_data.p'):
    """
Reads data from TPXO files and interpolates them onto the ROMS grid

Inputs:
------
  -lon, lat: ndarray type objects of (n,m) shape containing coordinates in deg.
  - harmonic_names: list of harmonic constituents, list
Options:
-------
  - pickle_name: name of the pickle file backup
  - tpxofile_data: dictionnary type object containing the paths of 
                        the tpxo netcdf data files.
    """

    lonmin1 = lon.min()
    lonmax1 = lon.max()
    latmin1 = lat.min()
    latmax1 = lat.max()

    if (lonmin1 < 0):
        lon = np.mod(lon, 360.0)

    if ((lonmin1 > 0) or (lonmax1 < 0)):
        pt1=(lon[0,0],lat[0,0])
        pt2=(lon[-1,0],lat[-1,0])
        pt3=(lon[-1,-1],lat[-1,-1])
        pt4=(lon[0,-1],lat[0,-1])
        bnd_vert = [pt1, pt2, pt3, pt4]
    else:
        pt1=(lon[0,0],lat[0,0])
        pt2=(lon[-1,0],lat[-1,0])
        pt3=(lon[-1,-1],lat[-1,-1])
        pt4=(lon[0,-1],lat[0,-1])
        pt5=(0.0,((lat[0,0]+lat[0,-1])/2.0))
        pt6=(0.0,((lat[-1,0]+lat[-1,-1])/2.0))
        pt7=(360.0,((lat[0,0]+lat[0,-1])/2.0))
        pt8=(360.0,((lat[-1,0]+lat[-1,-1])/2.0))
        bnd_vert_1= [pt1, pt2, pt7, pt8]
        bnd_vert_2= [pt6, pt6, pt3, pt4]

    lonmin = lon.min()
    lonmax = lon.max()
    latmin = lat.min()
    latmax = lat.max()

    grid30 = Dataset(tpxofile_data['grid'])

    X30 = grid30.variables['lon_z'][:]
    Y30 = grid30.variables['lat_z'][:]

    I30 = []
    J30 = []

    ##correction when grid spans over longitude 0
    if ((lonmin1 > 0) or (lonmax1 < 0)):
        #for i in range(0, len(X6)):
        #    if ((X6[i]<=(lonmax+0.5))and(X6[i]>=(lonmin-0.5))):
        #        I6.append(i)
        #for j in range(0, len(Y6)):
        #     if ((Y6[j]<=(latmax+0.5))and(Y6[j]>=(latmin-0.5))):
        #        J6.append(j)
        for i in range(0, len(X30)):
            if ((X30[i]<=(lonmax+0.5))and(X30[i]>=(lonmin-0.5))):
                I30.append(i)
        for j in range(0, len(Y30)):
             if ((Y30[j]<=(latmax+0.5))and(Y30[j]>=(latmin-0.5))):
                J30.append(j)
    else :
        #for i in range(0, len(X6)):
        #    if (((X6[i]<=(lonmax1 + 0.5))and(X6[i]>=0.0))\
        #    or ((X6[i]<=360.0)and(X6[i]>=np.mod(lonmin1-0.5,360.0)))):
        #        I6.append(i)
        #for j in range(0, len(Y6)): 
        #    if ((Y6[i,j]<=(latmax+0.5))and(Y6[i,j]>=(latmin-0.5))):
        #        J6.append(j)
        for i in range(0, len(X30)):
            if (((X30[i]<=(lonmax1 + 0.5))and(X30[i]>=0.0))\
            or ((X30[i]<=360.0)and(X30[i]>=np.mod(lonmin1-0.5,360.0)))):
                I30.append(i)
        for j in range(0, len(Y30)): 
            if ((Y30[i,j]<=(latmax+0.5))and(Y30[i,j]>=(latmin-0.5))):
                J30.append(j)            
 
    coord = ['z', 'u', 'v']
    var = ['h', 'u', 'v']
    TPXO = {}
    test=''

    for harmo in harmonic_names: 
        TPXO[harmo.lower()] = {}
        for v in var:
            TPXO[harmo.lower()][v] = {}

    for harmonic in harmonic_names:
        harmonic = harmonic.lower()
        print("Extracting {} components".format(harmonic))
        #print ('K=' + str(K))

        elev = Dataset(tpxofile_data['elev'][harmonic])
        vel = Dataset(tpxofile_data['vel'][harmonic])
        files = [elev,vel,vel]
        if harmonic in ['ms4', 'mn4', 'mf', 'mm']:
            grid = grid30
            I=I30
            J=J30
        else:
            grid = grid30
            I=I30
            J=J30
        
        for istep in range(0, len(var)):
            #print ('istep=' + str(istep))

            X1 = files[istep].variables['lon_' + coord[istep]][:]
            Y1 = files[istep].variables['lat_' + coord[istep]][:]
            [Y,X]=np.meshgrid(Y1,X1)
            H = grid.variables['h' + coord[istep]][:]
            Z = np.vectorize(complex)(files[istep].variables[var[istep] + 'Re'][:], \
                                      files[istep].variables[var[istep] + 'Im'][:])

            x = np.zeros((len(I), len(J)))
            y = np.zeros((len(I), len(J)))
            h = np.zeros((len(I), len(J)))
            z = np.zeros((len(I), len(J)))            
            z = np.vectorize(complex)(z, z)

            (row, col) = x.shape
            for i in range(0, row):
                for j in range(0, col):
                    if var[istep] == 'h':
                        x[i,j] = X[(I[i]),J[j]]
                        y[i,j] = Y[(I[i]),J[j]]                
                        h[i,j] = H[(I[i]),J[j]]
                        z[i,j] = Z[(I[i]),J[j]]/1000.0
                    else:                       
                        x[i,j] = X[I[i],J[j]]
                        y[i,j] = Y[I[i],J[j]]                
                        h[i,j] = H[I[i],J[j]]
                        z[i,j] = Z[I[i],J[j]]/10000.0

            x_pts, y_pts = x.flatten(), y.flatten()
            points = np.vstack((x_pts,y_pts)).T

            #correction when grid spans over longitude 0
            if ((lonmin1 > 0) or (lonmax1 < 0)):
                mask = Path(bnd_vert).contains_points(points)
                mask = mask.reshape((len(I), len(J)))
            else:
                mask_1 = Path(bnd_vert_1).contains_points(points)
                mask_1 = mask_1.reshape((len(I), len(J)))   
                mask_2 = Path(bnd_vert_2).contains_points(points)
                mask_2 = mask_2.reshape((len(I), len(J)))
                mask = (mask_1 + mask_2)
                mask = np.putmask(mask, mask>1, 1)             

            TPXO[harmonic][var[istep]] = {'x': x, 'y': y, 'h': h, 'mask': mask, 'z': z}

            if var[istep] != 'h':
                TPXO[harmonic][var[istep]]['z'][:,:] = z[:,:] / h[:,:]


    #Save TPXO file in pickle
    print("Pickling tidal data im {} ".format(pickle_name))
    pickle.dump(TPXO, open(pickle_name, "wb"))  
      
    return TPXO



###################################### Main function #######################################
def tpxo2roms(t0, ncROMSgridname, harmonic_names, tpxo_dict_paths,
              pickle_name='tpxo8_data.p', ncoutFilename='pcse_tides_3d_tpxo9_year.nc',ndays=397):
    """
Prepares a tidal forcing file for ROMS from TPXO tidal model (OSU) version 8.
The TPXO8  gathers 1/30th degree and 1/6th degree data.

Inputs:
------
  - t0: time for which phase is computed, "datetime" object
  - ncROMSgridname: name of the ROMS grid file, string
  - harmonic_names: list of harmonic constituents, list
Options:
-------
  - pickle_name: name of the back-up/restart file for tpxo8 info, string
  - ncoutFilename: name of the *.nc tide forcing file, string
  - ndays: Approximate anticipated length of model run (days), floats
    """
    #ROMS grid info
    roms_grid = Dataset(ncROMSgridname) 
    lon_rho = roms_grid.variables['lon_rho'][:]
    lat_rho = roms_grid.variables['lat_rho'][:]
    mask_rho = roms_grid.variables['mask_rho'][:]
    mask_psi = roms_grid.variables['mask_psi'][:]
    (M, L) = mask_rho.shape

    #Magic numbers for TPXO (Solar Doodson Numbers, Ref. Phase, Speed degrees/hour)
    tpxo_harmonics={}
    tpxo_harmonics['MM'] = [0, 1, 0, -1, 0, 0, 0, 0.5443747]
    tpxo_harmonics['MF'] = [0, 2, 0, 0, 0, 0, 0, 1.00980331]
    tpxo_harmonics['Q1'] = [1, -3, 1, 1, 0, 0, 270, 13.3986607]
    tpxo_harmonics['O1'] = [1, -2, 1, 0, 0, 0, 270, 13.9430351]
    tpxo_harmonics['P1'] = [1, 0, -1, 0, 0, 0, 270, 14.9589310]
    tpxo_harmonics['S1'] = [1, 0, 1, 0, 0, 0, 180, 15.0000000]
    tpxo_harmonics['K1'] = [1, 0, 1, 0, 0, 0, 90, 15.0410690]
    tpxo_harmonics['2N2'] = [2, -4, 2, 2, 0, 0, 0, 27.8953548]
    tpxo_harmonics['N2'] = [2, -3, 2, 1, 0, 0, 0, 28.4397297]
    tpxo_harmonics['M2'] = [2, -2, 2, 0, 0, 0, 0, 28.9841042]
    tpxo_harmonics['S2'] = [2, 0, 0, 0, 0, 0, 0, 30.0000000]
    tpxo_harmonics['K2'] = [2, 0, 2, 0, 0, 0, 0, 30.0821381]
    tpxo_harmonics['MN4'] = [4, -5, 4, 1, 0, 0, 0, 57.4238319]
    tpxo_harmonics['M4'] = [4 ,-4, 4, 0, 0, 0, 0, 57.9682083]
    tpxo_harmonics['MS4'] = [4, -2, 2, 0, 0, 0, 0, 58.9841042]

    #TPXO grid info
    tpxo=read_tpxo_data(lon_rho, lat_rho, harmonic_names, tpxo_dict_paths, pickle_name=pickle_name)
    
    #Magic numbers the way ROMS sees them
    roms_periods=[]
    vdeg={}
    for key in harmonic_names:
        roms_periods.append(360.0 / tpxo_harmonics[key][-1])
        vdeg[key] = vphase(t0, tpxo_harmonics[key])

    (F_FAC, U_FAC) = tpxo_nodal_factors(t0, ndays, harmonic_names)

    #Extract tide info from TPXO and put on rho grid
    N_HARMONICS = len(harmonic_names)
    ZAMP = np.zeros((N_HARMONICS, M, L))
    ZPHA = np.zeros((N_HARMONICS, M, L))
    UAMP = np.zeros((N_HARMONICS, M, L))
    UPHA = np.zeros((N_HARMONICS, M, L))
    VAMP = np.zeros((N_HARMONICS, M, L))
    VPHA = np.zeros((N_HARMONICS, M, L))
    MAJOR = np.zeros((N_HARMONICS, M, L))
    MINOR = np.zeros((N_HARMONICS, M, L))
    ECCENTRICITY = np.zeros((N_HARMONICS, M, L))
    INCLINATION = np.zeros((N_HARMONICS, M, L))
    PHASE = np.zeros((N_HARMONICS, M, L))
    
    for K in range(0, len(harmonic_names)):
        key = harmonic_names[K]
        print("Interpolating {} amplitudes".format(key))
        ZI = interp_tpxo(tpxo[key.lower()]['h'], lon_rho, lat_rho, mask_rho)
        ZAMP[K, :, :] = np.abs(ZI) * F_FAC[key] #element-wise multiplication
        ZPHA[K, :, :] = np.mod(-np.angle(ZI) * 180.0 / np.pi - U_FAC[key] - vdeg[key], 360.0)
        
        print("Interpolating {} amplitudes".format(key))
        UI = interp_tpxo(tpxo[key.lower()]['u'], lon_rho, lat_rho, mask_rho)
        UAMP[K, :, :] = np.abs(UI) * F_FAC[key] #element-wise multiplication
        UPHA[K, :, :] = np.mod(-np.angle(UI) * 180.0 / np.pi - U_FAC[key] - vdeg[key], 360.0)

        print("Interpolating {} amplitudes".format(key))
        VI = interp_tpxo(tpxo[key.lower()]['v'], lon_rho, lat_rho, mask_rho)
        VAMP[K, :, :] = np.abs(VI) * F_FAC[key] #element-wise multiplication
        VPHA[K, :, :] = np.mod(-np.angle(VI) * 180.0 / np.pi - U_FAC[key] - vdeg[key], 360.0)

        [MAJ, ECC, INC, PHA] = ampha2ellip(UAMP[K, :, :], UPHA[K, :, :], VAMP[K, :, :], VPHA[K, :, :])

        MAJOR[K, :, :] = MAJ
        ECC[np.isnan(ECC)] = 0.0
        MINOR[K, :, :] = MAJ * ECC
        #ECCENTRICITY[K, :, :] = ECC
        INCLINATION[K, :, :] = INC
        PHASE[K, :, :] = PHA
        
    #MINOR = MAJOR * ECCENTRICITY
    #MINOR = MAJ * ECC

    write_forcing(M, L, harmonic_names, roms_periods, ZAMP, ZPHA, UAMP, UPHA, VAMP, VPHA, MINOR, MAJOR, INCLINATION, PHASE, ncoutFilename)

###################################### Sub functions #######################################

def tpxo_nodal_factors(t0, dnum, Names):
    """
    Compute the f and u factors for the harmonics listed in cell array Names.
    Only those which are in the TPXO model are evaluated. They are evaluated at time dnum.

    See Table xxvi in A.T. Doodson (1928) 'On the Analysis of Tidal Observations'
    Philosophical Transactions of the Royal Society of London. Series A, Vol. 227
    J. Luick, Thanksgiving Day, 2011, Adelaide
        """

    end_time =datetime2datenum(t0) + dnum/2.0
    origin_time = datetime2datenum(datetime.datetime(1900, 1, 1))
    T = (end_time + 0.5 - origin_time) / 36525.0
    VN=np.mod(360.0*((0.719954-5.372617*T)+(0.000006*T*T)),360.0)
    if VN < 0.0: VN += 360
    VN = VN * np.pi / 180.0
    #VN = 3.3889

    #Coefficients
    CN = np.cos(VN)
    C2N = np.cos(2.0*VN)
    C3N = np.cos(3.0*VN)
    SN = np.sin(VN)
    S2N = np.sin(2.0*VN)    
    S3N = np.sin(3.0*VN)
    
    #Assign values for f and u of nonsolar constituents
    #(e.g. Doodson Table XXVI with u*pi/180)
    F = {}
    U = {}
    F['MM'] = 1.0 - (0.1300 * CN) + (0.0013 * C2N)
    U['MM'] = 0.0

    F['MF'] = 1.0429 + (0.4135 * CN) - (0.004 * C2N)
    U['MF'] = (-0.4143 * SN) + (0.0468 * S2N) - (0.0066 * S3N)
    
    F['O1'] = 1.0089 + (0.1871 * CN) - (0.0147 * C2N) + (0.0014 * C3N)
    U['O1'] = 0.1885 * SN - 0.0234 * S2N + 0.0033 * S3N

    F['K1'] = 1.0060 + (0.1150 * CN) - (0.0088 * C2N) + (0.0006 * C3N)
    U['K1'] = (-0.1546 * SN) + (0.0119 * S2N) - (0.0012 * S3N)

    F['M2'] = 1.0004 - (0.0373 * CN) + (0.0002 * C2N)
    U['M2'] = -0.0374 * SN

    F['K2'] = 1.0241 + (0.2863 * CN) + (0.0083 * C2N) - (0.0015 * C3N)
    U['K2'] = (-0.3096 * SN) + (0.0119 * S2N) - (0.0007 * S3N)

    #Dependent values
    F['Q1'] = F['O1']
    U['Q1'] = U['O1']
    F['N2'] = F['M2']
    U['N2'] = U['M2']
    F['2N2'] = F['M2']
    U['2N2'] = U['M2']
    F['MN4'] = F['M2'] ** 2.0
    U['MN4'] = 2.0 * U['M2']
    F['M4'] = F['M2'] ** 2.0
    U['M4'] = 2.0 * U['M2']
    F['MS4'] = F['M2']
    U['MS4'] = U['M2']

    #Assign F_Fac and U_Fac
    F_FAC = {}
    U_FAC = {}
    K = 0
    for name in Names:
        #if F.has_key(name):
        if name in F:
            F_FAC[name] = F[name]
            U_FAC[name] = np.mod(U[name] * 180.0 / np.pi, 360.0)
        else:
            #if f and u are not assigned, they are probably solar terms
            F_FAC[name] = 1.0
            U_FAC[name] = 0.0

    return (F_FAC, U_FAC)

def datetime2datenum(dt):
    """"
    Makes the conversion Python (e.g. datetime) and Matlab (e.g. datenum) time objects
    """    
    ord = dt.toordinal()
    mdn = dt + datetime.timedelta(days = 366)
    frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / \
           (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

def vphase(t0, dn_list):
    """
    Compute equilibrium phases in accordance with Cartwright "tidal analysis - a retrospect", 1982, pp170 - 188 in "Time series methods in hydrosciences,
    A.H.el-Shaarawi and S.R.Esterby (eds), Elsevier
    """
    end_time =datetime2datenum(t0)
    origin_time = datetime2datenum(datetime.datetime(1900, 1, 1))
    t = (end_time + 0.5 - origin_time) / 36525.0
    t_hour = np.mod(end_time, 1.0) * 24.0

    DN_NBR = dn_list[0:7]
    DN_PHA = dn_list[6]
    DN_SPD = dn_list[7]

    VS = np.mod(360.0 * (0.751206 + (1336.855231 * t) - (0.000003 * t * t)) , 360.0)
    VH = np.mod(360.0 * (0.776935 + ( 100.002136 * t) + (0.000001 * t * t)) , 360.0)
    VP = np.mod(360.0 * (0.928693 + (  11.302872 * t) - (0.000029 * t * t)) , 360.0)
    VN = np.mod(360.0 * (0.719954 - (   5.372617 * t) + (0.000006 * t * t)) , 360.0)
    VP1= np.mod(360.0 * (0.781169 + (   0.004775 * t) + (0.000001 * t * t)) , 360.0)

    for truc in [VS, VH, VP, VN, VP1]:
        if truc < 0: truc += 360.0

    Vdeg = (t_hour * DN_SPD) + (VS * DN_NBR[1]) + (VH * DN_NBR[2]) + \
           (VP * DN_NBR[3]) + (VN * DN_NBR[4]) + (VP1 * DN_NBR[5]) + DN_PHA

    return np.mod(Vdeg, 360.0)

def interp_tpxo(TPXOvar, lon, lat, mask):
    """
    Extract TPXO data and interpolate onto ROMS grid

    Inputs:
    VAR: hRE, hIm, uRE, uIm, vRe, or vIm
    lon, lat: 2D arrays of longitude and latitude to interpolate to
    mask: array of 1s and 0s corresponding to wet and dry points of lon & lat

    Output: VARinterp (VAR on lon, lat grid)

    Untested for domains extending across dateline - probably bogus there.
    TPXO files must be on the python path

    Created by J. Luick November 2011
    Modified by C James October 2013
    Transposed to python by T. Roc November 2013
    """

    iswet = (mask == 1).flatten()
    (I, J) = lon.shape

    lon = np.mod(lon, 360.0)

    x = TPXOvar['x'].flatten()
    y = TPXOvar['y'].flatten()
    z = TPXOvar['z'].flatten()
    h = TPXOvar['h']
    m = np.greater(np.logical_and(TPXOvar['mask'],h), 0).flatten()

    L = x.shape[0]
    orig = np.zeros((L, 2))
    orig[:, 0] = x
    orig[:, 1] = y   
 
    N = x[m].shape[0]
    known = np.zeros((N, 2))
    known[:, 0] = x[m]
    known[:, 1] = y[m]

    M = lon.flatten()[iswet].shape[0]
    asked = np.zeros((M, 2))
    asked[:, 0] = lon.flatten()[iswet]
    asked[:, 1] = lat.flatten()[iswet]   
    
    interpol = interpolate.NearestNDInterpolator(known, z[m])
    #interpol = interpolate.LinearNDInterpolator(known, z[m])
    #interpol = interpolate.CloughTocher2DInterpolator(known, z[m])
    z_interpol = interpol(orig)
    #interpol = interpolate.NearestNDInterpolator(orig, z_interpol)
    interpol = interpolate.LinearNDInterpolator(orig, z_interpol)
    #interpol = interpolate.CloughTocher2DInterpolator(orig, z_interpol)
    z_interpol = interpol(asked)

    VARinterp = np.zeros(iswet.shape).astype(complex)
    VARinterp[iswet] = z_interpol
    VARinterp = VARinterp.reshape((I, J))

    return VARinterp

def ampha2ellip(Au, PHIu, Av, PHIv):
    """
    Convert tidal amplitude and phase lag (ap-) parameters into tidal ellipse
    (e-) parameters. Please refer to ep2app for its inverse function.

    Usage:

    (SEMA,  ECC, INC, PHA, w)=app2ep(Au, PHIu, Av, PHIv, plot_demo)

    where:

        Au, PHIu, Av, PHIv are the amplitudes and phase lags (in degrees) of
        u- and v- tidal current components. They can be vectors or
        matrices or multidimensional arrays.
    http://www.acasports.co.uk/product_info.php?products_id=155
        plot_demo is an optional argument, when it is supplied as an array
        of indices, say [i j k l], the program will plot an  ellipse
        corresponding to Au(i,j, k, l), PHIu(i,j,k,l), Av(i,j,k,l), and
        PHIv(i,j,k,l);

        Any number of dimensions are allowed as long as your computer
        resource can handle.

        SEMA: Semi-major axes, or the maximum speed;
        ECC:  Eccentricity, the ratio of semi-minor axis over
            the semi-major axis; its negative value indicates that the ellipse
            is traversed in clockwise direction.
        INC:  Inclination, the angles (in degrees) between the semi-major
            axes and u-axis.
        PHA:  Phase angles, the time (in angles and in degrees) when the
            tidal currents reach their maximum speeds,  (i.e.
            PHA=omega*tmax).

            These four e-parameters will have the same dimensionality
            (i.e., vectors, or matrices) as the input ap-parameters.

        w:    Optional. If it is requested, it will be output as matrices
            whose rows allow for plotting ellipses and whose columns are
            for different ellipses corresponding columnwise to SEMA. For
            example, plot(real(w(1,:)), imag(w(1,:))) will let you see
            the first ellipse. You may need to use squeeze function when
            w is a more than two dimensional array. 
    """           
    
    #Assume the input phase lags are in degrees and convert them in radians.
    PHIu = PHIu / 180.0 * np.pi
    PHIv = PHIv / 180.0 * np.pi
    
    #Make complex amplitudes for u and v 
    u = Au * np.exp(-1j * PHIu)
    v = Av * np.exp(-1j * PHIv) 

    #Calculate complex radius, amplitudes and angles
    wp = (u + 1j * v) / 2.0        #for anticlockwise circle
    wm = np.conj(u - 1j * v) / 2.0 #for clockwise circles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)

    #Calculate e-parameters (ellipse parameters)
    SEMA = Wp + Wm              #Semi  Major Axis, or maximum speed
    SEMI = Wp - Wm              #Semin Minor Axis, or minimum speed
    ECC = np.divide(SEMI,SEMA)           #Eccentricity
    #ECC = SEMI / SEMA           #Eccentricity
    PHA = (THETAm - THETAp) / 2 #Phase angle
    INC = (THETAm + THETAp) / 2 #Inclination

    #Convert to degrees for output
    PHA = PHA / np.pi * 180.0
    INC = INC / np.pi * 180.0

    #Convention conversion for PHA, and INC
    id = PHA < 0
    PHA[id] = PHA[id] + 360.0
    id = INC < 0
    INC[id] = INC[id] + 360.0

    return SEMA, ECC, INC, PHA



def write_forcing(M, L, tpxo_harmonics, roms_periods, ZAMP, ZPHA, UAMP, UPHA, VAMP, VPHA, MINOR, MAJOR, INCLINATION, PHASE, ncoutFilename='output.nc'):
    """
Create and write in the netcdf forcing file
    """
    print("Creating {} netcdf file".format(ncoutFilename))
    #if os.path.isfile(ncoutFilename):
    #    CREATE_NEW = False
    #else:
    #	CREATE_NEW = True
    SIZE = len(tpxo_harmonics)
    nc = Dataset(ncoutFilename, 'w')
    nc.Description = 'ROMS forcing'
    nc.Author = 'Dr Thomas Roc'
    nc.Created = datetime.datetime.now().isoformat()
    nc.type = 'ROMS FRC file'
    nc.createDimension('xi_rho', L)
    nc.createDimension('eta_rho', M)
    nc.createDimension('tide_period', SIZE)

    nc.createVariable('tide_period', 'f8', 'tide_period')
    nc.variables['tide_period'].units = 'hours'
    nc.variables['tide_period'][:] = roms_periods

    nc.createVariable('tide_Eamp', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Eamp'].units = 'meters'
    nc.variables['tide_Eamp'][:] = ZAMP

    nc.createVariable('tide_Ephase', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Ephase'].units = 'degrees'
    nc.variables['tide_Ephase'][:] = ZPHA 

    nc.createVariable('tide_Uamp', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Uamp'].units = 'meter second-1'
    nc.variables['tide_Uamp'][:] = UAMP

    nc.createVariable('tide_Uphase', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Uphase'].units = 'degrees'
    nc.variables['tide_Uphase'][:] = UPHA
 
    nc.createVariable('tide_Vamp', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Vamp'].units = 'meter second-1'
    nc.variables['tide_Vamp'][:] = VAMP

    nc.createVariable('tide_Vphase', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Vphase'].units = 'degrees'
    nc.variables['tide_Vphase'][:] = VPHA

    nc.createVariable('tide_Cmin', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Cmin'].units = 'meter second-1'
    nc.variables['tide_Cmin'][:] = MINOR

    nc.createVariable('tide_Cmax', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Cmax'].units = 'meter second-1'
    nc.variables['tide_Cmax'][:] = MAJOR

    nc.createVariable('tide_Cangle', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Cangle'].units = 'degrees'
    nc.variables['tide_Cangle'][:] = INCLINATION

    nc.createVariable('tide_Cphase', 'f8', ('tide_period', 'eta_rho', 'xi_rho'))
    nc.variables['tide_Cphase'].units = 'degrees'
    nc.variables['tide_Cphase'][:] = PHASE


if __name__ == '__main__':
    
    #harmonic_names = ['K1','K2','M2','N2','O1','P1','Q1','S2','M4']
    harmonic_names = ['MM', 'MF', 'Q1', 'O1', 'P1', 'K1', 'S1', 'N2', 'M2', 'S2', 'K2', \
                      '2N2', 'MN4', 'M4', 'MS4']
    
    t0=datetime.datetime(2019,7,1,12,00,00)
    ndays = 550 # de 1/Jul/2019 -- 31/Dez/2020
    # 
    ncROMSgrdname = '/home/danilo/phd/thesis/chapter3_driving/projects/swa/deproas_v2002/sanagu2019/grid_SWA5km.nc'
    ncoutFilename = '/home/danilo/phd/thesis/chapter3_driving/projects/swa/deproas_v2002/sanagu2019/tides_SWA5km_SANAGU2019_2020.nc'
    tpxo2roms(t0, ncROMSgrdname, harmonic_names, ncoutFilename=ncoutFilename, ndays=ndays)