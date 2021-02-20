# ROMS TOOLS

## 1. Usage

### 1.1 Available functions
  This package uses pyroms generates:

 * regular grid (rotation option is available)
 * initial conditions
 * lateral boundary conditions
 * atmospheric forcing file
 * tpxo forcing file

The forcing files are not provided in this package, but are listed beow

### 1.2 Required datasets

Roms tools interpolates into ROMS netcdf input files information from :

* a. [Global Ocean Physics Reanalysis](https://resources.marine.copernicus.eu/?option=com_csw&task=results?option=com_csw&view=details&product_id=GLOBAL_REANALYSIS_PHY_001_030) (initial condition/boundary files)
* b. [ERA5 single levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form) (surface forcing)
* c. [ERA5 pressure levels](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form) (surface forcing)
* d. [TPXO](https://www.tpxo.net/) (tides)
* e. [GEBCO](https://www.gebco.net/data_and_products/gridded_bathymetry_data/) (bathymetry)
* f. Global Ocean Physics Reanalysis netcdf (GLO-MFC-001_030_mask_bathy.nc)

The file director is organized as (the letters indicate where you should add downloaded files):
```
.
├── data
│   ├── atmos_model <-- (b,c)
│   ├── aux_data    <-- (e,f)
│   ├── ocean_model <-- (a)
│   ├── roms_files  <-- files generated by this package should be stored here
│   ├── tpxo_files  <-- (d)
│   └── gridid.txt  <-- this is a txt file used to configure pyroms
└── scripts
    ├── era_glorys  
    └── tpxo

```

### 1.3 Running scripts

  After getting the datasets, configure the make_*.py according to your needs. You need to edit and export the path to data/gridid.txt (required by pyroms).
TO-DO: document this step


These are the main scripts we are using here.

 1. roms grid: make_regular_grid.py
 2. mask reference file: make_ref_mask.py
 3. weight files for interpolation: make_remaps_weights_file.py
 4. initial conditions: make_ic_files.py
 5. boundary conditions: make_bdry_files.py
 6. surface forcing: era2roms.py  <-- refactor era netcdfs to a  format roms can read


## 2. Installation

These tools are based on PyROMS, which is not simple to compile. For this reason we provide a Dockerfile. The following commands will allow you to run a docker container with pyroms and work  with pyroms_tools:

```
git clone https://github.com/CoastalHydrodynamicsLab/pyroms_tools.git
cd pyroms_tools

#building your container
sudo docker build -t pyroms_tools .

#after installation:
PREFIX=${PWD}
export UID=$(id -u)
sudo docker run -it --user=$UID -v $PREFIX/data:/home/lhico/data -v $PREFIX/scripts:/home/lhico/scripts pyroms

```

If you need to use the Graphic User Interface (GUI), you need to set up the paths from your machine to the docker container. In my case (Ubuntu 18.04), I did the following:

```
PREFIX=/your/pyroms_lhico/directory
sudo docker run -it \
   -v $PREFIX/data:/home/lhico/data \
   -v $PREFIX/scripts2:/home/lhico/scripts \
   -e DISPLAY=$DISPLAY \
   -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    pyroms
```

We don't include the installation of X11 in the Dockerfile, because it depends on hardware. A more in-depth explanation may be found [here](https://stackoverflow.com/questions/25281992/alternatives-to-ssh-x11-forwarding-for-docker-containers)
