modelsvarbs = {
    'glorys':{
        'global_reanalysis_phy_001_031':{
            'lon': 'longitude',
            'lat': 'latitude',
            'mlt': 'mlotst_%s',  # mixed layer depth
            'temp':'thetao_%s',  # potential temperature
            'sal': 'so_%s',  # practical salinity
            'u':   'uo_%s',  # eastward velocity
            'v':   'vo_%s',  # northward velocity
            'ssh': 'zos_%s',  # sea surface height
        }
    }
}


area_params = {
    'wSAtl0': {
        'lonW': -70,
        'lonE':   0,
        'latS': -60,
        'latN':   0,
        'xrange': (200, 500),
        'yrange': (300, 500)
    },
    'pcse0': {
        'lonW': -55,
        'lonE': -35,
        'latS': -35,
        'latN':  -20
    },
    'pcse1': {
        'lonW': -55,
        'lonE': -30,
        'latS': -35,
        'latN':  -15,
        'xrange': (10, 291),
        'yrange': (10, 231)
    },
    'Florence': {
        'lonW': -85,
        'lonE': -47,
        'latS': 18,
        'latN':  40
    }
}
