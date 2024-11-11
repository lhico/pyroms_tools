from pyroms_tools.utils import edit_mask as em
from pyroms_tools.utils import utils as ut
import xarray as xr
from scipy.ndimage import label
import numpy as np


def largest_contiguous_ocean_mask(matrix):
    # Create a mask where True represents non-NaN elements
    mask = matrix == 1
    
    # Label contiguous blocks of True values
    labeled_array, num_features = label(mask)
        
        # Count the size of each contiguous region
    # Label contiguous regions
    labeled_mask, num_features = label(mask)

    # Find the largest contiguous region
    largest_region_label = np.argmax(np.bincount(labeled_mask.flat)[1:]) + 1  # ignore the 0 label (background)

    # Create a new matrix that only includes the largest region
    largest_region = np.where(labeled_mask == largest_region_label, matrix, False)

    
    return largest_region

if __name__ == '__main__':
    # ds = xr.open_dataset('/mnt/34c919f6-6617-49de-a20f-05e4186230b9/Dropbox/trabalho_irado/Northeastern/other/roms_grid01_smooth.nc')
    ds = xr.open_dataset('roms_grid01_smooth.nc')
    ds.mask_rho.values = largest_contiguous_ocean_mask(ds.mask_rho.values)

    em.main(ds.mask_rho.values, ds.lat_rho.values, ds.lon_rho.values)

    ds = ut.update_mask(ds)