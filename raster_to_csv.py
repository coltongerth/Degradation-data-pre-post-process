import os
import rasterio
import asyncio
import warnings
# import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio.mask
from shapely.geometry import box
# from shapely.geometry import shape
from scipy import stats as st
from datetime import datetime
from rasterio.mask import mask
from dask.delayed import delayed
# from dotenv import load_doten
from IPython.display import display
from rasterio.windows import Window
from rasterio.profiles import DefaultGTiffProfile
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED, ALL_COMPLETED

def main():
        
    # IF YOU HAVE MEMORY CRASHING ISSUES LOWER THIS NUMBER (KEEP IN POWERS OF 2: 2, 4, 8, 16, 32 ...)
    BLOCKSIZE = 256

    # INPUTS HERE
    zone_name = "BpsMskR2Fin"
    boundary_name = "BpsMskR2Fin"

    # DONT WORRY ABOUT THESE
    zone_raster_path = f"{zone_name}.tif"
    post_cut_zone_raster_path = f"./data/{zone_name}/{zone_name}_{boundary_name}_cut.tif"
    boundary_raster_path = f"{boundary_name}.tif"
    data_raster_path = f"./data/{zone_name}/{boundary_name}_rpms_stack.tif"
    dummy_path = "./test.tif"
    data_type = "float32"

    def raster_to_csv(raster, file_path):
        """Convert a raster dataset to CSV, ensuring valid values are extracted."""
        with rasterio.open(raster) as src:
            # Read metadata and print for debugging
            print("Metadata:", src.meta)

            # Read the first band
            band1 = src.read(1)

            # Print minimum and maximum values for debugging
            print("Min value:", band1.min())
            print("Max value:", band1.max())

            # Check the NoData value
            nodata = src.nodata
            print("NoData value:", nodata)

            # Read rows and cols indices
            rows, cols = np.indices(band1.shape)

            # Transform the indices to coordinates
            x_coords, y_coords = src.transform * (cols, rows)

            # Flatten arrays
            x_coords = x_coords.flatten()
            y_coords = y_coords.flatten()
            values = band1.flatten()

            # Filter out NoData values if NoData is not None
            if nodata is not None:
                mask = values != nodata
                x_coords = x_coords[mask]
                y_coords = y_coords[mask]
                values = values[mask]

            # Print a sample of values for debugging
            print("Sample values:", values[:10])

            # Create a DataFrame and save to CSV
            df = pd.DataFrame({
                'X': x_coords,
                'Y': y_coords,
                'Zone': values
            })

            # Save to CSV
            df.to_csv(file_path, index=False)
            print(f"CSV saved to {file_path}")

    def write_multiband_raster_to_single_csv(raster, file_name):
        """Utility function to write multiband raster data to a single CSV file."""
        with rasterio.open(raster) as src:
            all_bands_data = {}

            # Get the coordinates once since they are the same for all bands
            rows, cols = np.indices((src.height, src.width))
            x_coords, y_coords = src.transform * (cols, rows)
            x_coords = x_coords.flatten()
            y_coords = y_coords.flatten()

            # Store the coordinates in the DataFrame dictionary
            all_bands_data['X'] = x_coords
            all_bands_data['Y'] = y_coords

            # Iterate over each band
            for band_index in range(1, src.count + 1):  # src.count gives the number of bands
                band = src.read(band_index)  # Read each band
                values = band.flatten()

                # Add the band's values to the dictionary with a key for each band
                all_bands_data[f'value_band_{band_index}'] = values

            # Convert the dictionary to a DataFrame
            df = pd.DataFrame(all_bands_data)

            # drop rows where all bands are 0
            df = df[(df != 0).any(axis=1)]

            # Save the DataFrame to a CSV file
            df.to_csv(file_name, index=False)

    # async def raster_stacker(in_ds, out_ds, bounds):
    def raster_stacker(id, in_ds, out_ds, bounds):

        with rasterio.open(in_ds) as src_ds:
            win = src_ds.window(
                bottom=bounds.bottom,
                right=bounds.right,
                top=bounds.top,
                left=bounds.left,
            )
            print(f"in: {in_ds} || {win}")
            input_data = src_ds.read(1, window=win)
            out_ds.write_band(id, input_data)



    # with rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS", verbose=2, GOOGLE_APPLICATION_CREDENTIALS=os.getenv("GOOGLE_APPLICATION_CREDENTIALS", path_to_credentials)):
    with rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS", verbose=2):

        boundary_ds = rasterio.open(boundary_raster_path)
        bounds = boundary_ds.bounds
        # bounds = [-114.174042,42.633959,-112.858429,43.608239]

        profile = boundary_ds.profile
        profile.update(
            blockxsize=BLOCKSIZE,
            blockysize=BLOCKSIZE,
            tiled=True,
            compress="DEFLATE",
            # predictor=2,
            BIGTIFF="Yes",
            dtype=data_type

        )

        od = f"./data/{zone_name}"
        if not os.path.exists(od):
            os.makedirs(od)

        if os.path.exists(post_cut_zone_raster_path):
            print(f"cut zone data {post_cut_zone_raster_path} exists.")
        else:
            # Open the larger raster and cut it using the bounds of the smaller raster
            with rasterio.open(zone_raster_path) as dx:
                # Calculate the window from the smaller raster's bounds
                win = dx.window(left=bounds.left, bottom=bounds.bottom,
                                right=bounds.right, top=bounds.top)

                # Read the data from the larger raster within the window
                dat = dx.read(window=win)
                # Adjust the profile for the output file
                profile = dx.profile
                profile.update({
                    'height': dat.shape[1],
                    'width': dat.shape[2],
                    'transform': rasterio.windows.transform(win, dx.transform)
                })

                # Write the cut raster to a new file
                with rasterio.open(post_cut_zone_raster_path, 'w', **profile) as dst:
                    dst.write(dat)


        raster_to_csv(post_cut_zone_raster_path, f'{zone_name}_{boundary_name}_cut.csv')


        files = [f"https://storage.googleapis.com/fuelcast-public/rpms/{y}/rpms_{y}.tif" for y in range(1984, 2024) if y != 2012]

        profile.update(count=len(files))

        print("Stacking raster")

        stack_path = data_raster_path

        if os.path.exists(stack_path):
            print(f"Stacked raster {stack_path} already exists.")
        else:
            with rasterio.open(stack_path, "w", **profile) as dst:
                print(f"out: {dst} || {dst.bounds}")

                for id, layer in enumerate(files, start=1):
                    print(f"in: {layer}")
                    raster_stacker(id, layer, dst, bounds)
                # raster_stacker(1, files[0], dst, bounds)


        write_multiband_raster_to_single_csv(stack_path, f'{boundary_name}_rpms_stack_cut.csv')

        zones_df = pd.read_csv(f'{zone_name}_{boundary_name}_cut.csv')
        zones_df = zones_df[zones_df['value'] != -32768]
        zones_df.rename(columns={'x': 'X', 'y': 'Y'}, inplace=True)
        zones_df.to_csv(f'{zone_name}_{boundary_name}_cut.csv', index=False)

        new_cut_rpms_filtered = pd.read_csv(f'{boundary_name}_rpms_stack_cut.csv')

        def replace_zeros_with_row_mean(df):
            # Extract only the 'band' columns for processing
            band_columns = [col for col in df.columns if col.startswith('band')]

            # Iterate over each row
            for index, row in df.iterrows():
                # Calculate the mean of non-zero values in the row
                mean_value = row[band_columns][row[band_columns] != 0].mean()

                # Replace zeroes with the calculated mean
                df.loc[index, band_columns] = row[band_columns].replace(0, mean_value)

            return df

        new_cut_rpms_filtered_mean_filled = replace_zeros_with_row_mean(new_cut_rpms_filtered)

        new_cut_rpms_filtered_mean_filled.to_csv(f'{boundary_name}_rpms_stack_cut_filtered_mean_filled.csv', index=False)




if __name__ == "__main__":
    main()
