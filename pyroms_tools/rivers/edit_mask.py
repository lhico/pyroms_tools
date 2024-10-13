import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Create a large numpy boolean matrix (e.g., representing lat/lon grid)
matrix = np.zeros((180, 720), dtype=bool)

# Define latitude and longitude ranges corresponding to the matrix
lat = np.linspace(90, -90, matrix.shape[0])  # Latitudes from 90 to -90
lon = np.linspace(-180, 180, matrix.shape[1])  # Longitudes from -180 to 180

# Editing mode toggle
editing_mode = False  # Initially, editing mode is off

# Function to find the nearest matrix index for lat/lon
def find_nearest_idx(array, value):
    return np.abs(array - value).argmin()

# Function to toggle matrix values on click
def on_click(event):
    global editing_mode, pc, fig
    if editing_mode and event.inaxes:  # Only change values if editing mode is active
        # Convert click coordinates (lon, lat) to matrix indices
        lon_clicked, lat_clicked = event.xdata, event.ydata
        x_idx = find_nearest_idx(lon, lon_clicked)
        y_idx = find_nearest_idx(lat, lat_clicked)

        # Toggle the value at the clicked position
        matrix[y_idx, x_idx] = not matrix[y_idx, x_idx]

        # Update the pcolormesh object with the new data
        pc.set_array(matrix.ravel())
        pc.set_clim([matrix.min(), matrix.max()])  # Adjust color limits if needed
        plt.draw()

# Function to toggle editing mode on/off when a key is pressed
def on_key(event):
    global editing_mode
    if event.key == 'e':  # Press 'e' to toggle editing mode
        editing_mode = not editing_mode
        print(f"Editing mode {'ON' if editing_mode else 'OFF'}")

# Set up the figure with Cartopy's projection
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_global()  # Show the whole world map

# Add initial land contours and coastlines
ax.add_feature(cfeature.LAND, edgecolor='black')
ax.coastlines()

# Overlay the matrix on top of the Cartopy map using pcolormesh
pc = ax.pcolormesh(lon, lat, matrix, cmap='gray', transform=ccrs.PlateCarree(), shading='auto', alpha=0.5, edgecolors='black', linewidth=0.5,)

ax.grid(True, which='both', color='black', linestyle='-', linewidth=0.5)  # Draw gridlines

# Connect the click event to the on_click function
fig.canvas.mpl_connect('button_press_event', on_click)

# Connect the key press event to the on_key function for toggling editing mode
fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()
