from utils import click_map as cp
# from utils import utils as ut
import xarray as xr
import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
import sys

ffile = sys.argv[1] #if len(sys.argv) == 2 else "sbb_grid_roms.nc"

nc = xr.open_dataset(ffile)
plt.close('all')
tracker = cp.IndexTracker(nc, ttype='mask')
while tracker.wait:
    fig, ax = tracker.get_figure()
    # fig.tight_layout()
    plt.pause(2)
    tracker.show()
    plt.pause(1)

tracker.dataset.to_netcdf('modified_grid.nc')