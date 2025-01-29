# Degradation-data-pre-post-process
The python script version of the pre-processing and post processing of raster data to use with the R Degradation Script

How to Run

Step 1: Download the required python libraries. In your commandline, ensure you are in the directory with the raster_to_csv.py script and run

python3 -m pip install -r requirements.txt

Step 2: Ensure that the zones raster is in the same local directory it should look as follows in your file explorer

/forRScript
bpsgt100ku.tif
raster_to_csv.py
README.md
requirements.txt
result_csv_to_raster.py

Step 3: Run the python script to form an RPMS raster stack within the bounds of our zones raster:

python3 raster_to_csv.py

This should result in two csv files in the forRScript folder, you'll need to take these and place them in the local directory of where you run the R Script for Degradation. If you run the R Script from within this folder, all you need to do is add the folder location "./forRScript" to the front of the two input csvs.

Step 4: Run the R Script Deg_with_zoning.R. I run this in Rstudio and just highlight all the lines of code in the GUI and click the run button. This should result in a csv when all is said and done called: 

comparison_analysis_zones_rpms_full_with_good_zone_res.csv

Step 5: Run result_csv_to_raster.py:

python3 result_csv_to_raster.py

This should produce four rasters (t_mean.tif, t_slope.tif, p_mean.tif, p_slope.tif), if that's the case then mission success!!!

