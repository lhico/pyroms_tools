micromamba install -c conda-forge python=3.11

micromamba install -c conda-forge dask netCDF4 ipython --yes
micromamba install -c conda-forge esmpy xarray numpy shapely cf_xarray sparse numba bottleneck cartopy --yes
micromamba install -c conda-forge xesmf --yes
micromamba install -c conda-forge scikit-fmm pyproj tqdm gsw scikit-build-core cmake ninja scikit-build lpsolve55 --yes

# pyroms installation
git clone https://github.com/lhico/pyroms

pip install build

cd pyroms/pyroms
pip install .

cd ../pyroms_toolbox
pip install .

cd ../bathy_smoother
pip install .


python -m build
pip install .

cd vertical_interpolation
pip install .
cd ..

cd extrapolate
pip install .
cd ..
