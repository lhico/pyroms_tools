# Preprocessing roms files

## 1. Introduction
This set of scripts were prepared to create roms input files. These scripts are based on the pyroms and xesmf modules and covers:
* orthogonal rectangular rotatated grid  creation;
* initial condition creation;
* boundary conditions creation;
* tide forcing;

There are three files that configures the grid and interpolation details. **Ideally, for the preprocessing, we should only change the information within these configuration files instead of working directly with the python scripts**. The configuration files are:

```
.
├── configs 
│   ├── grid_config_esmf.txt
│   ├── grid_config_pyroms.txt
│   ├── gridid.txt
(...)
```

Within these configuration files we have paths, grid configurations, variable names maps and other information that used by the scripts.

### **1.1 Installation (pyroms)**

These tools are based on PyROMS, which is not simple to compile. For this reason we provide a docker image and a Dockerfile. There are two ways to install it. The following commands will allow you to run a docker container with pyroms and work  with pyroms_tools:

1) Pulling from a docker image (recommended)
```
git clone https://github.com/lhico/pyroms_tools.git
cd pyroms_tools

#pulling your container
docker pull oteldks/pyroms_tools:latest

#run container after installation:
PREFIX=${PWD}
export UID=$(id -u)
sudo docker run -t -it --user=$UID -v $PREFIX:/home/lhico/pyroms_tools  oteldks/pyroms_tools:latest
```

2) Building from a docker file (not recommended).
```
git clone https://github.com/lhico/pyroms_tools.git
cd pyroms_tools

#building your container
sudo docker build -t pyroms_tools .

#run container after installation:
PREFIX=${PWD}
export UID=$(id -u)
sudo docker run -it --user=$UID -v $PREFIX:/home/lhico/pyroms_tools   pyroms_tools

```

If you have problems with the command `python somescript.py` within the container, please type `source ~/.bashrc` to activate conda environment.

<!-- During the building step, the following error might appear:

```bash
sudo: effective uid is not 0, is /usr/bin/sudo on a file system with the 'nosuid' option set or an NFS file system without root privileges?
```

If this is your case, then you must use the ```Dockerfile.v2``` instead. To do this, try the following sequence of commands:

```bash
#building your container
sudo docker build -t pyroms_tools . -f Dockerfile.v2

#after installation:
PREFIX=${PWD}
export UID=$(id -u)
sudo docker run -it --user=$UID -v $PREFIX:/home/lhico/pyroms_tools   pyroms_tools

``` -->

### **1.2 GUI interface with Docker**

If you need to use the Graphic User Interface (GUI), you need to set up the paths from your machine to the docker container. In my case (Ubuntu 18.04), I did the following:

```
PREFIX=/your/pyroms_lhico/directory
sudo docker run -it  \
    -v $PREFIX:/home/lhico/pyroms_tools \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    oteldks/pyroms_tools: latest
```

We don't include the installation of X11 in the Dockerfile, because it depends on hardware. A more in-depth explanation may be found [here](https://stackoverflow.com/questions/25281992/alternatives-to-ssh-x11-forwarding-for-docker-containers)


### **1.3 Installation xesmf

In order to have a working xesmf environment using miniconda, use the following commands.

```
conda create -n xesmf_teste python=3.7
conda activate xesmf_teste
conda install -c conda-forge xesmf
conda install -c conda-forge xesmf esmpy=8.0.0
conda install -c conda-forge dask netCDF4
conda install -c conda-forge matplotlib cartopy
```

## 2. Scripts structure

The main scripts are placed at the `scripts` and `scripts_xesmf` directory as follows:

```
.
├── scripts
│   ├── CGrid_TPXO8   # directory
│   ├── utils         # directory
│   ├── clear
│   ├── make_bdry_files.py
│   ├── make_bdry_remap_weight_files.py
│   ├── make_grid_prototype.py
│   ├── make_grid_.py
│   ├── make_grid_smooth.py
│   ├── make_ic_file.py
│   ├── make_ic_file_replace.py
│   ├── make_ic_remap_weight_files.py
│   ├── make_tpxo_remap_weight_file.py
│   ├── make_tpxo_tide.py
│   └── README.md
└── scripts_xesmf
    ├── make_forcing_files.py
    ├── make_ic_file_replace.py
    ├── nested_adjust_grids_bathymetry.py
    ├── nested_replace_finegrid_bathy.py
    └── utils -> ../scripts/utils/  # symlink
```

Our scripts are used to prepare the following files:

1. grids
2. initial conditions
3. boundary files
4. tides
5. nesting adjustments

Notice that these scripts assume the dataset used as a reference for the interpolation is in an Arakawa A-grid (https://en.wikipedia.org/wiki/Arakawa_grids). Their paths are set within the configuration files.


The files in `scripts` should be used within a pyroms environments - we provide a docker image with a pyroms environment, since this package may be rather difficult to install.

The files in `scripts_xesmf` require (do not use pyroms docker image):
- xarray
- numpy
- scipy
- xesmf
- cartopy

If you are using the nested scripts, please run the following line in `scripts/utis`:

```
f2py --verbose -c -m extrapolate fill.f90
```

If you are using nudging scripts, please run the following line in `scripts/utils`:
```
f2py -c -m interpolation interpolation.f90
```


### **2.1 Grids**

We need to create roms grid:

a. create a prototype grid: `make_grid_prototype.py` (pyroms, optional)

b. create the grid: `make_grid_.py` (pyroms)

c. smooth the grid: `make_grid_smooth.py` (pyroms, optional).
Grid smoothing should adjust the values rx0 (0\<rx0\<3) and rx1 (3\<rx0\<8) to avoid hydrostatic inconsistencies due to the s-coordinate system. These values are printed when roms simulation is started. (optional) 

Everything is linked by the dictionary configuration file `grid_config_pyroms.txt`. Notice that item `c.` is just an usefule example and there are matlab scripts also used for smoothing. You must run the grid with idealized stratified conditions to check if 'motionless currents' are too strong. They appear due to the numerical inconsistensies that arise from sigma-coordinates and should not present values over a few centimeters per second.

### **2.2 Initial conditions**

After creating the grid file, you'll need to interpolate the initial conditions onto it. We have three steps in this case:
    
a. create empty ic file : `make_ic_file.py` (pyroms)

b. interpolate the information onto empty ic_file: `make_ic_file_replace.py` (xesmf)

### **2.3 Boundary files**

After creating the grid file, you'll need to interpolate the boundary conditions onto boundary condition files. We have two steps in this case:

a. create interpolation weights: `make_bdry_files_remap_weight_file.py` (pyroms)

b. interpolate the information onto roms grid: `make_bdry_file.py` (pyroms)


### **2.4 Tide files**


After creating the grid file, you'll need to interpolate the tpxo data onto the grid. The tides are not forced by the boundary, nor the surface. It is a field to avoid sponges that may be set close to the boundary points. We have two steps in this case:

a. create interpolation weights: `make_tpxo_files_remap_weight_file.py` (pyroms)

b. interpolate the information onto roms grid: `make_tpxo_file.py` (pyroms)


**ATTENTION!!!!!**: tpxo files are huge and we don't want to copy them repeatedly. Docker containers don't play well with symbolic links, so the workaround to have tpxo files in the right place is to map the volume directly from your computer onto a container directory. Below is an example:

```
TPXOPATH=/home/otel/Dropbox/trabalho_irado/2021/postdoc/2021_data/roms_files/tpxo
PREFIX=/path/to/pyroms_tools
 sudo docker run -it -v $TPXOPATH:/home/lhico/pyroms_tools/data/raw/tpxo -v $PREFIX:/home/lhico/pyroms_tools pyroms_tools
```


### Interpolation

Since we use pyroms to interpolate information onto roms grid, it is necessary to create interpolation weights - `*_remap_weight_files.py` - prior to the interpolation files.



<!---
## Nesting procedure (not completed)

1 - Create coarse grid  CRS
2 - Create refined grid REC1 (coarse2fine_group.m part A matlab)
3 - Create refined grid SRC1 (make_grid_refine.py): this will be used as a source of information. REC1 domain must be contained in SRC1 domain, but doesn't need to be the same asa CRS
4 - Adjust bathymetries and masks to grids (REC1,2) (make_grid_refine*_interpolation.py)
5 - create remap files
6 - create initial conditions files
7 - reinterpolate initial condition files
8 - run the final part in c
9 - take the grid files and create the contact file (coarse2fine_group.m part B matlab)
10 - if you are nesting you need to adjust the bathymetry between the fields. it is possible to make the adjustment by adapting the script script_aux/adjust_bathy. Remeber, the contour lines must be overlapping on the boundaries, otherwise the model will generate weird gradients in the contact points
-->
