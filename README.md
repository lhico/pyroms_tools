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
├── data
│   ├── grid_config_esmf.txt
│   ├── grid_config_pyroms.txt
│   ├── gridid.txt
│   └── roms_files
(...)
```

Within these configuration files we have paths, grid configurations, variable names maps and other information that used by the scripts.

## 2. Scripts structure

The main scripts are placed at the `scripts` directory as follows:

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
│   ├── make_grid_refine1_interpolation.py
│   ├── make_grid_refine2_interpolation.py
│   ├── make_grid_refine.py
│   ├── make_grid_smooth.py
│   ├── make_ic_file.py
│   ├── make_ic_file_replace.py
│   ├── make_ic_remap_weight_files.py
│   ├── make_tpxo_remap_weight_file.py
│   ├── make_tpxo_tide.py
│   └── README.md
│   
```

Our scripts are used to prepare the following files:

1. grids
2. initial conditions
3. boundary files
4. tides

Notice that these scripts assume the dataset used as a reference for the interpolation is in an Arakawa A-grid (https://en.wikipedia.org/wiki/Arakawa_grids). Their paths are set within the configuration files.

### **2.1 Grids**

We need to create roms grid:

a. create a prototype grid: `make_grid_prototype.py` (pyroms, optional)

b. create the grid: `make_grid_py` (pyroms)

c. smooth the grid: `make_grid_smooth.py` (pyroms, optional).
Grid smoothing should adjust the values rx0 (0\<rx0\<3) and rx1 (3\<rx0\<8) to avoid hydrostatic inconsistencies due to the s-coordinate system. These values are printed when roms simulation is started. (optional) 

Everything is linked by the dictionary configuration file `grid_config_pyroms.txt`. Notice that item `c.` is just an usefule example and there are matlab scripts also used for smoothing. You must run the grid with idealized stratified conditions to check if 'motionless currents' are too strong. They appear due to the numerical inconsistensies that arise from sigma-coordinates and should not present values over a few centimeters per second.

### **2.2 Initial conditions**

After creating the grid file, you'll need to interpolate the initial conditions onto it. We have three steps in this case:
    
a. create interpolation weights: `make_ic_files_remap_weight_file.py` (pyroms)

b. interpolate the information onto roms grid: `make_ic_file.py` (pyroms)

c. reinterpolate the information onto roms grid: `make_ic_file_replace.py` (xesmf)

The third step is necessary in these scripts because the interpolation using pyroms created spurious horizontal TS gradients in idealized cases where the TS fields were  horizontally homogenous. For this reason I rewrote an interpolation script with xesmf that corrected the issue. **Warning: close to the coast I needed to interpolate the information with a nearest-neighbor approach, which is NOT general. This quick-fix worked in my case, but it needs further work.** 

### **2.3 Boundary files**

After creating the grid file, you'll need to interpolate the boundary conditions onto boundary condition files. We have two steps in this case:

a. create interpolation weights: `make_bdry_files_remap_weight_file.py` (pyroms)

b. interpolate the information onto roms grid: `make_bdry_file.py` (pyroms)


### **2.4 Tide files**

After creating the grid file, you'll need to interpolate the tpxo data onto the grid. The tides are not forced by the boundary, nor the surface. It is a field to avoid sponges that may be set close to the boundary points. We have two steps in this case:

a. create interpolation weights: `make_tpxo_files_remap_weight_file.py` (pyroms)

b. interpolate the information onto roms grid: `make_tpxo_file.py` (pyroms)



### Interpolation

Since we use pyroms to interpolate information onto roms grid, it is necessary to create interpolation weights - `*_remap_weight_files.py` - prior to the interpolation files.




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
