import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
from rasterio.transform import from_origin
from google.colab import drive

data = pd.read_csv('comparison_analysis_zones_rpms_full_with_good_zone_res.csv')

# CHANGE THIS IF YOUR CSV HAS X,Y or LATITUDE,LONGITUDE
lat_long = False

if lat_long == False:
    x_col = "X"
    y_col = "Y"
else:
    x_col = "longitude"
    y_col = "latitude"



# want t test, slope test, t pvalue, slope pvalue
bands = {
    "t(mean)":"t_mean",
    "t(slope)":"t_slope",
    "p(|t|<T)(mean)":"p_mean",
    "p(|t|<T)(slope)":"p_slope"
}

pixel_size = 0.0002694945852358564

# Define the extent of the data
x_min, x_max = data[x_col].min(), data[x_col].max()
y_min, y_max = data[y_col].min(), data[y_col].max()

# Create grid resolution
x_res = int((x_max - x_min) / pixel_size) + 1
y_res = int((y_max - y_min) / pixel_size) + 1

for band in bands.keys():
    # Initialize an empty array for the raster (NaNs for no data)
    raster_data = np.full((y_res, x_res), np.nan)

    # Map DataFrame coordinates to the grid
    for index, row in data.iterrows():
        # Use round() and clip to avoid out-of-bounds errors
        col = np.clip(round((row[x_col] - x_min) / pixel_size), 0, x_res - 1)
        row_ = np.clip(round((y_max - row[y_col]) / pixel_size), 0, y_res - 1)  # Y-axis is inverted in rasters

        # Set the value in the raster data (e.g., p-value for slope)
        raster_data[row_, col] = row[band]

    # Define transformation (mapping pixel coordinates to geographic coordinates)
    transform = from_origin(x_min, y_max, pixel_size, pixel_size)

    with rasterio.open(
        f'{bands[band]}_raster.tif', 'w',
        driver='GTiff',
        height=raster_data.shape[0],
        width=raster_data.shape[1],
        count=1,
        dtype=raster_data.dtype,
        crs='EPSG:4326',
        transform=transform
    ) as dst:
        dst.write(raster_data, 1)