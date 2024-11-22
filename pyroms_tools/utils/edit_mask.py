import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def find_nearest_idx(lon_array, lat_array, lon_value, lat_value):
    """
    Find the indices of the nearest point in the 2D lon and lat arrays to the given lon and lat values.
    
    Parameters:
    - lon_array (np.ndarray): 2D array of longitudes.
    - lat_array (np.ndarray): 2D array of latitudes.
    - lon_value (float): Reference longitude.
    - lat_value (float): Reference latitude.
    
    Returns:
    - (int, int): Indices (j, i) of the nearest point.
    """
    # Calculate the squared distance to avoid computing the square root
    dist_squared = (lon_array - lon_value)**2 + (lat_array - lat_value)**2
    j, i = np.unravel_index(dist_squared.argmin(), dist_squared.shape)
    return j, i

def on_click(event, matrix, lat, lon, pc):
    if event.inaxes:  # Only change values if editing mode is active
        # Convert click coordinates (lon, lat) to matrix indices
        lon_clicked, lat_clicked = event.xdata, event.ydata
        y_idx, x_idx = find_nearest_idx(lon, lat, lon_clicked, lat_clicked)

        # Toggle the value at the clicked position
        matrix[y_idx, x_idx] = not matrix[y_idx, x_idx]

        # Update the pcolormesh object with the new data
        pc.set_array(matrix.ravel())
        pc.set_clim([matrix.min(), matrix.max()])  # Adjust color limits if needed
        plt.draw()

def on_key(event, editing_mode):
    if event.key == 'e':  # Press 'e' to toggle editing mode
        editing_mode[0] = not editing_mode[0]
        print(f"Editing mode {'ON' if editing_mode[0] else 'OFF'}")

def main(matrix, lat, lon):
    editing_mode = [False]  # Initially, editing mode is off

    # Set up the figure with Cartopy's projection
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_global()  # Show the whole world map

    # Add initial land contours and coastlines
    ax.add_feature(cfeature.LAND, edgecolor='black')
    ax.coastlines()

    # Overlay the matrix on top of the Cartopy map using pcolormesh
    pc = ax.pcolormesh(lon, lat, matrix, cmap='gray', transform=ccrs.PlateCarree(), shading='auto', alpha=0.5, edgecolors='black', linewidth=0.5)

    ax.grid(True, which='both', color='black', linestyle='-', linewidth=0.5)  # Draw gridlines

    # Connect the click event to the on_click function
    fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, matrix, lat, lon, pc))

    # Connect the key press event to the on_key function for toggling editing mode
    fig.canvas.mpl_connect('key_press_event', lambda event: on_key(event, editing_mode))

    plt.show()
    return matrix

if __name__ == "__main__":
    # Create a large numpy boolean matrix (e.g., representing lat/lon grid)
    matrix = np.zeros((180, 720), dtype=bool)

    # Define latitude and longitude ranges corresponding to the matrix
    lat = np.linspace(90, -90, matrix.shape[0])  # Latitudes from 90 to -90
    lon = np.linspace(-180, 180, matrix.shape[1])  # Longitudes from -180 to 180
    xm,ym = np.meshgrid(lon,lat)    

    main(matrix, xm, ym)
    