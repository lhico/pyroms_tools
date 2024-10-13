micromamba install -c conda-forge python=3.11

micromamba install -c conda-forge dask netCDF4 ipython
micromamba install -c conda-forge esmpy xarray numpy shapely cf_xarray sparse numba bottleneck
micromamba install -c conda-forge xesmf
micromamba install -c conda-forge scikit-fmm pyproj tqdm

# pyroms installation
git clone https://github.com/lhico/pyroms

micromamba install -c conda-forge scikit-build-core cmake ninja
micromamba install scikit-build
pip install build
micromamba install -c conda-forge lpsolve55

cd pyroms/pyroms
pip install .

cd ../pyroms_toolbox
pip install .

cd ../bathy_smoother
pip install .


python -m build
pip install .

