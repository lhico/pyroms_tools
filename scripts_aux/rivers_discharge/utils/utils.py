import xarray as xr 
import re

def river_netcdf(riverid, s_rho, time,
                 xpos, epos, sign, direction,
                 Vshape, transport, temp, salt,
                 lon, lat):
    """Generates the netcdf river input for ROMS

    Args:
        xpos (array like):
            river XI-position
            dims=('river')
        epos (array like):
            river ETA-position
            dims=('river')
        sign (array like):
            controls the direction of the flow acrross a u-face or v-face source
            positive value is source, negative value is a sink
            dims=('river')
        direction (array like):
            When using option LuvSrc = T, the variable river_direction defines
            whether the river enters a cell through the u-face (river_direction = 0)
            or v-face (river_direction = 1) of the receiving cell.
            dims=('river')
        Vshape (array like): 
            Variable river_Vshape sets the fractional distribution of the
            river_transport among the vertical cells and must sum to 1 over
            the vertical. If the sum over k is not 1 then the actual volume flux will depart from the value given in river_transport.
            dims=('s_rho', 'river')
        transport (array like):
            river runoff vertically integrated mass transport
            dims=('river_time', 'river')
        temp (array like):
            river runoff potential temperature
            dims=('river_time', 's_rho', 'river')
        salt (array like):
            river runoff practical salinity
            dims=('river_time', 's_rho', 'river')

    Returns:
        xarray Dataset
    """

    river = {}
    river['long_name'] = "river runoff identification number" ;

    river_time = {}
    river_time['long_name'] = "river runoff time" ;
    # river_time['units'] = "days since 2001-01-01 00:00:00" ;

    river_direction = {}
    river_direction['long_name'] = "river runoff direction" ;
    river_direction['flag_values'] = "0, 1" ;
    river_direction['flag_meanings'] = "flow across u-face, flow across v-face" ;
    river_direction['LwSrc_True'] = "flag not used" ;

    river_Xposition = {}
    river_Xposition['long_name'] = "river XI-position" ;
    river_Xposition['LuvSrc_meaning'] = "i point index of U or V face source/sink" ;
    river_Xposition['LwSrc_meaning'] = "i point index of RHO center source/sink" ; 

    river_Eposition = {}
    river_Eposition['long_name'] = "river ETA-position" ;
    river_Eposition['LuvSrc_True_meaning'] = "j point index of U or V face source/sink" ;
    river_Eposition['LwSrc_True_meaning'] = "j point index of RHO center source/sink" ;

    river_transport = {}
    river_transport['long_name'] = "river runoff vertically integrated mass transport" ;
    river_transport['units'] = "meter3 second-1" ;
    river_transport['positive'] = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell" ;
    river_transport['negative'] = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell" ;
    river_transport['time'] = "river_time" ;

    river_Vshape    = {}
    river_Vshape['long_name'] = "river runoff mass transport vertical profile" ;
    river_Vshape['requires'] = "must sum to 1 over s_rho" ;

    river_temp      = {}
    river_temp['long_name'] = "river runoff potential temperature" ;
    river_temp['units'] = "Celsius" ;
    river_temp['time'] = "river_time" ;

    river_salt      = {}
    river_salt['long_name'] = "river runoff salinity" ;
    river_salt['time'] = "river_time" ;

    river_sign      = {}



    data_vars = {
        'river_Xposition': (('river'), xpos, river_Xposition),
        'river_Eposition': (('river'), epos, river_Eposition),
        'river_sign'     : (('river'), sign, river_sign),
        'river_direction': (('river'), direction, river_direction),
        'river_Vshape'   : (('s_rho', 'river'), Vshape, river_Vshape),
        'river_transport': (('river_time', 'river'), transport, river_transport),
        'river_temp'   : (('river_time', 's_rho', 'river'), temp, river_temp),
        'river_salt'   : (('river_time', 's_rho', 'river'), salt, river_salt),
        'lon': (('river'), lon),
        'lat': (('river'), lat),
    }

    coords = {
        'river_time': time,
        'river': riverid,
        's_rho': s_rho    }

    attrs = {
        'type': 'ROMS RIVERS file',
        'title': 'Testing',
        'source': 'Testing'
    }
    nc = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)
    return nc


def add_to_lists(pairs, i, j, sign, dir):
    x1, y1 = pairs[0]

    for it in range(1,len(pairs)):
        x2, y2 = pairs[it]

        if x2 > x1:
        # negative v-velocity
            i.append(x1)
            j.append(y1)
            sign.append(-1)
            dir.append(1)
        elif x1 > x2:
        # positive v-velocity
            i.append(x2)
            j.append(y1)
            sign.append(1)
            dir.append(1)
        elif y2 > y1:
        # positive u-velocity
            i.append(x1)
            j.append(y1)
            sign.append(1)
            dir.append(0)
        elif y1 > y2:
        # negative u-velocity
            i.append(x1)
            j.append(y2)
            sign.append(-1)
            dir.append(0)
        x1 = x2
        y1 = y2


def maskedge_points(landmask_boundary='maskedge.out'):
    # We need to parse the output of the maskedge program for two
    # different purposes:
    #  1. Create the rivers file for ROMS, at least the locations part.
    #  2. Create a scrip grid file for the river locations.
    # This routine will only do #1 (so far).

    # Read the landmask boundaries
    f = open('maskedge.out', 'r')
    pairs = []
    # Eat first line so we don't trigger the add_to_lists routine
    f.readline()

    # These are for the ROMS sources file
    i = []
    j = []
    sign = []
    dir = []

    #pdb.set_trace()

    for line in f:
        a, b, c = re.split('\s+', line)
        if a=='-10':
            # wrap up object
            add_to_lists(pairs, i, j, sign, dir)
        elif (a=='-1' or a=='-3'):
            # wrap up object
            add_to_lists(pairs, i, j, sign, dir)
            # start new object
            pairs = []
        else:
            pairs.append([int(a),int(b)])
            

    return pairs