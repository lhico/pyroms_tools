####################
# structure based on ROMS.jl (Alexander Barth), which was based on d_ecmwf2roms.m 
# from Arango and Wilkins
dt = 1 # in hour because it will be converted to seconds

variables_list = {
    'msl' : {
        'ECMWFlongname': 'mean_sea_level_pressure',
        'Vname': 'Pair',
        'accumulation': False,
        'outputName': 'Pair',
        'scale': 0.01, # Pa to mb
        'units': 'mb',
        'time': 'pair_time'
    },
    'u10' : {
        'ECMWFlongname': '10m_u_component_of_wind',
        'Vname': 'Uwind',
        'accumulation': False,
        'outputName': 'Uwind',
        'scale': 1.0,
        'units': 'm s-1',
        'time': 'wind_time'
    },
    'v10' : {
        'ECMWFlongname': '10m_v_component_of_wind',
        'Vname': 'Vwind',
        'accumulation': False,
        'outputName': 'Vwind',
        'scale': 1.0,
        'units': 'm s-1',
        'time': 'wind_time'
    },
    't2m' : {
        'ECMWFlongname': '2m_temperature',
        'Vname': 'Tair',
        'accumulation': False,
        'outputName': 'Tair',
        'scale': 1.0,
        'units': 'degC',
        'time': 'tair_time'
    },
    'tcc' : {
        'ECMWFlongname': 'total_cloud_cover',
        'Vname': 'cloud',
        'accumulation': False,
        'outputName': 'cloud',
        'scale': 1.0,
        'units': '',
        'time': 'cloud_time'
    },
    'msshf' : {
        'ECMWFlongname': 'surface_sensible_heat_flux',
        'Vname': 'sensible',
        'accumulation': True,
        'outputName': 'sensible',
        'scale': -1.0/(dt*3600),
        'units': 'W m-2',
        'time': 'sen_time'
    },
    'mslhf' : {
        'ECMWFlongname': 'surface_latent_heat_flux',
        'Vname': 'latent',
        'accumulation': True,
        'outputName': 'latent',
        'scale': -1.0/(dt*3600),
        'units': 'W m-2',
        'time': 'lhf_time'
    },
    'msnlwrf' : { # str
        'ECMWFlongname': 'surface_net_thermal_radiation',
        'Vname': 'lwrad',
        'accumulation': True,
        'outputName': 'lwrad',
        'scale': 1.0/(dt*3600),
        'units': 'W m-2',
        'time': 'wind_time'
    },
    'msdwlwrf' : { # strd
        'ECMWFlongname': 'surface_thermal_radiation_downwards',
        'Vname': 'lwrad_down',
        'accumulation': True,
        'outputName': 'lwrad_down',
        'scale': 1.0/(dt*3600),
        'units': 'W m-2'
    },
    'ssr' : {
        'ECMWFlongname': 'surface_net_solar_radiation',
        'Vname': 'swrad',
        'accumulation': True,
        'outputName': 'swrad',
        'scale': 1.0/(dt*3600),
        'units': 'W m-2'
    },
    'msnswrf' : {
        'ECMWFlongname': 'mean_surface_net_short_wave_radiation_flux',
        'Vname': 'swrad',
        'accumulation': True,
        'outputName': 'swrad',
        'scale': 1.0/(dt*3600),
        'units': 'W m-2'
    },
    'metss' : { # wess
        'ECMWFlongname': 'eastward_turbulent_surface_stress',
        'Vname': 'sustr',
        'accumulation': True,
        'outputName': 'sustr',
        'scale': 1.0/(dt*3600),
        'units': 'N m-2'
    },
    'mntss' : { # nsss
        'ECMWFlongname': 'northward_turbulent_surface_stress',
        'Vname': 'svstr',
        'accumulation': True,
        'outputName': 'svstr',
        'scale': 1.0/(dt*3600),
        'units': 'N m-2'
    },
    'mtpr' : {
        'ECMWFlongname': 'total_precipitation',
        'Vname': 'rain',
        'accumulation': True,
        'outputName': 'rain',
        'scale': 1.0/(dt*3600), #1000.0/(dt*3600)
        #'units': 'kg m-2 s-1'
    },
    'mer' : {
        'ECMWFlongname': 'evaporation',
        'Vname': 'e',
        'accumulation': True,
        'outputName': 'e',
        'scale': 1.0/(dt*3600), #1000.0/(dt*3600)
        #'units': 'kg m-2 s-1'
    },
    'sst' : {
        'ECMWFlongname': 'sea_surface_temperature',
        'Vname': 'sst',
        'accumulation': False,
        'outputName': 'sst',
        'scale': 1,
        'units': 'degC'
    },

    # variables to be derived from other variables
    'q' : {
        'ECMWFlongname': '2m_dewpoint_temperature',
        'Vname': 'Qair',
        'accumulation': False,
        'outputName': 'Qair',
        'scale': 1.0,
        'units': 'percentage'
    },
    'shflux' : {
        'ECMWFlongname': '',
        'Vname': 'shflux',
        'accumulation': True,
        'outputName': 'shflux',
        'scale': 1.0/(dt*3600),
        'units': 'W m-2'
    },
    'swflux' : {
        'ECMWFlongname': '',
        'Vname': 'swflux',
        'accumulation': True,
        'outputName': 'swflux',
        'scale': 1/1000, #-100. / (dt*3600.)*(24*3600), # (Arango's scripts)
        'units': 'm s-1'
    },
    'dQdSST' : {
        'ECMWFlongname': '',
        'Vname': 'dQdSST',
        'accumulation': True,
        'outputName': 'dQdSST',
        'scale': 1,
        'units': 'W m-2 degC-1'
    },
}